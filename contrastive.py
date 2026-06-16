import torch
import torch.nn as nn
import torch.optim as optim
# from store.code.make_datasets_ import *
# from models import *
import random
import numpy as np
import torch.nn.functional as F
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score, normalized_mutual_info_score, homogeneity_score
import itertools


class SupConLoss(nn.Module):
    """
    Supervised Contrastive Loss
    同类拉近，不同类拉远
    """

    def __init__(self, temperature=0.07):
        super().__init__()
        self.temperature = temperature

    def forward(self, features, labels):
        """
        features: [B, D]  (embedding)
        labels:   [B]
        """
        device = features.device
        B = features.size(0)

        # 1. L2 normalize
        # features = F.normalize(features, dim=1)

        # 2. 相似度矩阵 [B, B]
        sim = torch.matmul(features, features.T) / self.temperature

        # 3. mask
        labels = labels.contiguous().view(-1, 1)
        mask_pos = torch.eq(labels, labels.T).float().to(device)
        mask_self = torch.eye(B, device=device)
        mask_pos = mask_pos - mask_self  # 去掉自己

        # 4. log-softmax（分母）
        sim_max, _ = torch.max(sim, dim=1, keepdim=True)
        sim = sim - sim_max.detach()  # 数值稳定

        exp_sim = torch.exp(sim) * (1 - mask_self)
        log_prob = sim - torch.log(exp_sim.sum(dim=1, keepdim=True) + 1e-12)

        # 5. 对正样本求平均
        mean_log_prob_pos = (mask_pos * log_prob).sum(dim=1) / (
            mask_pos.sum(dim=1) + 1e-12
        )

        loss = -mean_log_prob_pos.mean()
        return loss


def train_encoder_classification_expand(E_rna, E_atac, Classifier, optimizer_E_rna, optimizer_E_atac, optimizer_fc, scheduler_E_rna, scheduler_E_atac, scheduler_fc, omics1_train_loader, omics1_test_id_loader, omics2_train_loader, omics2_test_id_loader, filtered_omics1_id_pseudo_loader, filtered_omics2_id_pseudo_loader, filtered_omics1_pseudo_ood_loader, num_epochs):
    E_rna.train()
    E_atac.train()
    Classifier.train()

    Loss_SupCon = SupConLoss(temperature=0.07)  # 0.07

    for epoch in range(num_epochs):
        running_loss = 0.0

        # 选长的那个作为主循环
        if len(omics1_train_loader) >= len(omics2_train_loader):
            main_loader = omics1_train_loader
            aux_iter = itertools.cycle(omics2_train_loader)
            main_is_rna = True
        else:
            main_loader = omics2_train_loader
            aux_iter = itertools.cycle(omics1_train_loader)
            main_is_rna = False

        # for (_, x_rna, y_pseudo_rna), (_, x_atac, y_pseudo_atac) in zip(omics1_train_loader, omics2_train_loader):
        for batch_main, batch_aux in zip(main_loader, aux_iter):
            if main_is_rna:
                (_, x_rna, y_pseudo_rna), (_, x_atac, y_pseudo_atac) = batch_main, batch_aux
            else:
                (_, x_atac, y_pseudo_atac), (_, x_rna, y_pseudo_rna) = batch_main, batch_aux

            x_rna = x_rna.view(x_rna.size(0), -1).cuda()  # 将图片展平为一维向量
            x_atac = x_atac.view(x_atac.size(0), -1).cuda()  # 将图片展平为一维向量

            y_pseudo_rna = y_pseudo_rna.view(-1)
            y_pseudo_atac = y_pseudo_atac.view(-1)

            y = torch.cat((y_pseudo_rna, y_pseudo_atac), dim=0).cuda()

            z_rna, _, _ = E_rna(x_rna)
            z_atac, _, _ = E_atac(x_atac)

            z = torch.cat((z_rna, z_atac), dim=0)

            logits = Classifier(z)

            loss_ce = F.cross_entropy(logits, y)
            loss_sc = Loss_SupCon(z, y)
            loss = loss_ce + loss_sc

            # 清空梯度
            optimizer_E_rna.zero_grad()
            optimizer_E_atac.zero_grad()
            optimizer_fc.zero_grad()

            loss.backward()  # 反向传播

            # 更新权重
            optimizer_E_rna.step()
            optimizer_E_atac.step()
            optimizer_fc.step()

            # 更新学习率
            scheduler_E_rna.step()
            scheduler_E_atac.step()
            scheduler_fc.step()

        # ACC_Train = evaluate_model(E_rna, E_atac, Classifier, omics1_train_loader, omics2_train_loader)

        ACC_Test = evaluate_model(E_rna, E_atac, Classifier, filtered_omics1_pseudo_ood_loader, filtered_omics2_id_pseudo_loader)

        # if ACC_Train > 99.0:
        #     break

        print(f'Epoch [{epoch+1}/{num_epochs}], Loss: {running_loss/len(omics1_train_loader):.4f}, LR: {scheduler_E_rna.get_last_lr()[0]:.6f}, Acc_test: {ACC_Test:.2f}')
        evaluate_model_with_for_novel_class(E_rna, Classifier, filtered_omics1_pseudo_ood_loader)  # filtered_omics1_pseudo_ood_loader

    return E_rna, E_atac, Classifier


def get_num_classes_from_loader(loader):
    labels = []
    for _, _, y in loader:
        y = y.view(-1)
        labels.append(y)
    labels = torch.cat(labels)

    return int(labels.max().item() + 1)


def expand_classifier(classifier, num_classes_new):
    """
    扩展线性分类头到 num_classes_new
    保留旧类别参数
    """
    assert hasattr(classifier, "cls"), "Classifier must have attribute `cls`"
    assert isinstance(classifier.cls, nn.Linear)

    old_linear = classifier.cls
    in_dim = old_linear.in_features
    num_classes_old = old_linear.out_features

    if num_classes_new <= num_classes_old:
        print("No expansion needed.")
        return classifier

    device = classifier.cls.weight.device
    dtype = classifier.cls.weight.dtype

    # 1. 新建更大的 classifier
    new_linear = nn.Linear(in_dim, num_classes_new).to(device=device, dtype=dtype)

    # 2. 拷贝旧权重
    with torch.no_grad():
        new_linear.weight[:num_classes_old] = old_linear.weight
        new_linear.bias[:num_classes_old] = old_linear.bias

        # # 3. 初始化新增类别（可选策略）
        # nn.init.normal_(
        #     new_linear.weight[num_classes_old:], mean=0.0, std=0.01
        # )
        # nn.init.constant_(
        #     new_linear.bias[num_classes_old:], 0.0
        # )

    classifier.cls = new_linear

    return classifier


def evaluate_model_with_for_novel_class(Encoder, Classifier, loader):
    Encoder.eval()
    Classifier.eval()

    true_label = []
    predicted_label = []
    with torch.no_grad():  # 关闭梯度计算
        for (y_true, x, _) in loader:
            # x_rna = x_rna.view(x_rna.size(0), -1).cuda()  # 将图片展平为一维向量
            # x_atac = x_atac.view(x_atac.size(0), -1).cuda()  # 将图片展平为一维向量
            x = x.cuda()

            y_true = y_true.view(-1)

            z, _, _ = Encoder(x)

            logits = Classifier(z)

            _, predicted = torch.max(logits.data, 1)

            true_label.append(y_true)
            predicted_label.append(predicted.cpu())

    true_label = torch.cat(true_label)
    predicted_label = torch.cat(predicted_label)

    ari = adjusted_rand_score(true_label, predicted_label)
    ami = adjusted_mutual_info_score(true_label, predicted_label)
    nmi = normalized_mutual_info_score(true_label, predicted_label)
    hom = homogeneity_score(true_label, predicted_label)

    print(f'Performance of Model on the Novel Class: ARI:{ari:.2f}, AMI:{ami:.2f}, NMI:{nmi:.2f}, HOM:{hom:.2f}')

