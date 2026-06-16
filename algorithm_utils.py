import torch
import random
import numpy as np
import torch.nn.functional as F
from scDiscovery_architecture import reparameterize
from data_split import BuildDataset
from contrastive import SupConLoss
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt


seed = 0
random.seed(0)
torch.manual_seed(seed)
np.random.seed(seed)
torch.cuda.manual_seed(seed)
rng = np.random.default_rng(seed)
torch.cuda.manual_seed_all(seed)


def train_encoder_classification(E_rna, E_atac, Classifier, optimizer_E_rna, optimizer_E_atac, optimizer_fc, scheduler_E_rna, scheduler_E_atac, scheduler_fc, omics1_train_loader, omics1_test_id_loader, omics2_train_loader, omics2_test_id_loader, num_epochs):
    E_rna.train()
    E_atac.train()
    Classifier.train()
    for epoch in range(num_epochs):
        running_loss = 0.0
        for (_, x_rna, y_rna), (_, x_atac, y_atac) in zip(omics1_train_loader, omics2_train_loader):
            # x_rna = x_rna.view(x_rna.size(0), -1).cuda()  # 将图片展平为一维向量
            # x_atac = x_atac.view(x_atac.size(0), -1).cuda()  # 将图片展平为一维向量
            x_rna = x_rna.cuda()
            x_atac = x_atac.cuda()

            y = torch.cat((y_rna, y_atac), dim=0).cuda()

            z_rna, _, _ = E_rna(x_rna)
            z_atac, _, _ = E_atac(x_atac)

            z = torch.cat((z_rna, z_atac), dim=0)

            logits = Classifier(z)

            loss = F.cross_entropy(logits, y)

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

        ACC_Train = evaluate_model(E_rna, E_atac, Classifier, omics1_train_loader, omics2_train_loader)

        ACC_Test = evaluate_model(E_rna, E_atac, Classifier, omics1_test_id_loader, omics2_test_id_loader)

        # if ACC_Train > 99.0:
        #     break

        print(f'Epoch [{epoch+1}/{num_epochs}], Loss: {running_loss/len(omics1_train_loader):.4f}, LR: {scheduler_E_rna.get_last_lr()[0]:.6f}, Acc_train: {ACC_Train:.2f}, Acc_test: {ACC_Test:.2f}')

    return E_rna, E_atac, Classifier


def train_encoder_classification_with_single_omics(E_rna, Classifier, optimizer_E_rna, optimizer_fc, scheduler_E_rna, scheduler_fc, omics1_train_loader, omics1_test_id_loader, num_epochs):
    E_rna.train()
    Classifier.train()
    for epoch in range(num_epochs):
        running_loss = 0.0
        for _, x_rna, y_rna in omics1_train_loader:
            x_rna = x_rna.cuda()
            y = y_rna.cuda()

            z_rna, _, _ = E_rna(x_rna)

            logits = Classifier(z_rna)

            loss = F.cross_entropy(logits, y)

            # 清空梯度
            optimizer_E_rna.zero_grad()
            optimizer_fc.zero_grad()

            loss.backward()  # 反向传播

            # 更新权重
            optimizer_E_rna.step()
            optimizer_fc.step()

            # 更新学习率
            scheduler_E_rna.step()
            scheduler_fc.step()

        ACC_Train = evaluate_model_with_single_omics(E_rna, Classifier, omics1_train_loader)

        ACC_Test = evaluate_model_with_single_omics(E_rna, Classifier, omics1_test_id_loader)

        # if ACC_Train > 99.0:
        #     break

        print(f'Epoch [{epoch+1}/{num_epochs}], Loss: {running_loss/len(omics1_train_loader):.4f}, LR: {scheduler_E_rna.get_last_lr()[0]:.6f}, Acc_train: {ACC_Train:.2f}, Acc_test: {ACC_Test:.2f}')

    return E_rna, Classifier


# 评估模型
def evaluate_model(E_rna, E_atac, Classifier, omics1_train_loader, omics2_train_loader, device):
    E_rna.eval()
    E_atac.eval()
    Classifier.eval()

    correct = 0
    total = 0
    with torch.no_grad():  # 关闭梯度计算
        for (_, x_rna, y_rna), (_, x_atac, y_atac) in zip(omics1_train_loader, omics2_train_loader):
            # x_rna = x_rna.view(x_rna.size(0), -1).cuda()  # 将图片展平为一维向量
            # x_atac = x_atac.view(x_atac.size(0), -1).cuda()  # 将图片展平为一维向量
            x_rna = x_rna.to(device)
            x_atac = x_atac.to(device)

            y = torch.cat((y_rna, y_atac), dim=0).to(device)

            z_rna, _, _ = E_rna(x_rna)
            z_atac, _, _ = E_atac(x_atac)

            z = torch.cat((z_rna, z_atac), dim=0)

            logits = Classifier(z)

            _, predicted = torch.max(logits.data, 1)

            total += logits.size(0)
            correct += (predicted == y).sum().item()
    accuracy = 100 * correct / total
    # print(f'Accuracy of the model on the test images: {accuracy:.2f}%')
    return accuracy


def evaluate_model_with_single_omics(E_rna, Classifier, omics1_loader, device):
    E_rna.eval()
    Classifier.eval()

    correct = 0
    total = 0
    with torch.no_grad():  # 关闭梯度计算
        for _, x_rna, y_rna in omics1_loader:
            x_rna = x_rna.to(device)
            y_rna = y_rna.to(device)

            z_rna, _, _ = E_rna(x_rna)

            logits = Classifier(z_rna)

            _, predicted = torch.max(logits.data, 1)

            total += logits.size(0)
            correct += (predicted == y_rna).sum().item()
    accuracy = 100 * correct / total
    # print(f'Accuracy of the model on the test images: {accuracy:.2f}%')
    return accuracy


def train_decoder_classification_single_omics_base(E_rna, Classifier, optimizer_E_rna, optimizer_fc, scheduler_E_rna, scheduler_fc, omics1_train_loader, omics1_test_id_loader, num_epochs, device):
    E_rna.train()
    Classifier.train()

    Loss_SupCon = SupConLoss(temperature=0.07)  # 0.07
    for epoch in range(num_epochs):
        total_loss_cls = 0.0
        total_rna_recon_loss = 0.0
        total_rna_kl_loss = 0.0
        for _, x_rna, y_rna in omics1_train_loader:
            x_rna = x_rna.to(device)

            y = y_rna.to(device)

            z_rna, mu_rna, logvar_rna = E_rna(x_rna)

            logits = Classifier(z_rna)

            loss_cls = F.cross_entropy(logits, y)

            # ---------------- decoder -----------------

            loss = 1.0*loss_cls  # + rna_kl_loss + rna_recon_loss + 0.5*loss_sc  #

            # 清空梯度
            optimizer_E_rna.zero_grad()
            optimizer_fc.zero_grad()

            loss.backward()  # 反向传播

            # 更新权重
            optimizer_E_rna.step()
            optimizer_fc.step()

            # 更新学习率
            scheduler_E_rna.step()
            scheduler_fc.step()

            total_loss_cls += loss_cls

        ACC_Train = evaluate_model_with_single_omics(E_rna, Classifier, omics1_train_loader, device)

        if omics1_test_id_loader is not None:
            ACC_Test = evaluate_model_with_single_omics(E_rna, Classifier, omics1_test_id_loader, device)
        else:
            ACC_Test = ACC_Train

        print(f'Epoch [{epoch+1}/{num_epochs}], CLS Loss: {total_loss_cls/len(omics1_train_loader):.4f}, RNA Recon Loss: {total_rna_recon_loss/len(omics1_train_loader):.4f}, RNA KL Loss: {total_rna_kl_loss/len(omics1_train_loader):.4f}, Acc_train: {ACC_Train:.2f}, Acc_test: {ACC_Test:.2f}')

        # if ACC_Train > 90.0:
        #     break

    return E_rna, Classifier


def train_decoder_classification(E_rna, E_atac, D_rna, D_atac, Classifier, optimizer_E_rna, optimizer_E_atac, optimizer_D_rna, optimizer_D_atac, optimizer_fc, scheduler_E_rna, scheduler_E_atac, scheduler_D_rna, scheduler_D_atac, scheduler_fc, omics1_train_loader, omics1_test_id_loader, omics2_train_loader, omics2_test_id_loader, num_epochs, device):
    E_rna.train()
    E_atac.train()
    D_rna.train()
    D_atac.train()
    Classifier.train()

    Loss_SupCon = SupConLoss(temperature=0.07)  # 0.07
    for epoch in range(num_epochs):
        total_loss_cls = 0.0
        total_rna_recon_loss = 0.0
        total_atac_recon_loss = 0.0
        total_rna_kl_loss = 0.0
        total_atac_kl_loss = 0.0
        for (_, x_rna, y_rna), (_, x_atac, y_atac) in zip(omics1_train_loader, omics2_train_loader):
            # x_rna = x_rna.view(x_rna.size(0), -1).cuda()  # 将图片展平为一维向量
            # x_atac = x_atac.view(x_atac.size(0), -1).cuda()  # 将图片展平为一维向量
            x_rna = x_rna.to(device)
            x_atac = x_atac.to(device)

            y = torch.cat((y_rna, y_atac), dim=0).to(device)

            z_rna, mu_rna, logvar_rna = E_rna(x_rna)
            z_atac, mu_atac, logvar_atac = E_atac(x_atac)

            re_z_rna = reparameterize(mu_rna, logvar_rna)
            re_z_atac = reparameterize(mu_atac, logvar_atac)

            re_z = torch.cat((re_z_rna, re_z_atac), dim=0).to(device)

            z = torch.cat((z_rna, z_atac), dim=0)

            logits = Classifier(z)

            loss_cls = F.cross_entropy(logits, y)

            # ---------------- decoder -----------------
            x_hat_rna = D_rna(re_z_rna)
            x_hat_atac = D_atac(re_z_atac)

            rna_recon_loss = F.mse_loss(x_hat_rna, x_rna, reduction="mean")
            atac_recon_loss = F.mse_loss(x_hat_atac, x_atac, reduction="mean")

            rna_kl_loss = -0.5 * torch.mean(1 + logvar_rna - mu_rna.pow(2) - logvar_rna.exp())
            atac_kl_loss = -0.5 * torch.mean(1 + logvar_atac - mu_atac.pow(2) - logvar_atac.exp())

            loss_sc = Loss_SupCon(F.normalize(z), y) + Loss_SupCon(F.normalize(re_z), y)

            loss = 0.1*loss_cls + rna_recon_loss + atac_recon_loss + rna_kl_loss + atac_kl_loss + 0.5*loss_sc

            # 清空梯度
            optimizer_E_rna.zero_grad()
            optimizer_E_atac.zero_grad()
            optimizer_D_rna.zero_grad()
            optimizer_D_atac.zero_grad()
            optimizer_fc.zero_grad()

            loss.backward()  # 反向传播

            # 更新权重
            optimizer_E_rna.step()
            optimizer_E_atac.step()
            optimizer_D_rna.step()
            optimizer_D_atac.step()
            optimizer_fc.step()

            # 更新学习率
            scheduler_E_rna.step()
            scheduler_E_atac.step()
            scheduler_D_rna.step()
            scheduler_D_atac.step()
            scheduler_fc.step()

            total_loss_cls += loss_cls
            total_rna_recon_loss += rna_recon_loss
            total_atac_recon_loss += atac_recon_loss
            total_rna_kl_loss += rna_kl_loss
            total_atac_kl_loss += atac_kl_loss

        ACC_Train = evaluate_model(E_rna, E_atac, Classifier, omics1_train_loader, omics2_train_loader, device)

        ACC_Test = evaluate_model(E_rna, E_atac, Classifier, omics1_test_id_loader, omics2_test_id_loader, device)

        print(f'Epoch [{epoch+1}/{num_epochs}], CLS Loss: {total_loss_cls/len(omics1_train_loader):.4f}, RNA Recon Loss: {total_rna_recon_loss/len(omics1_train_loader):.4f}, ATAC Recon Loss: {total_atac_recon_loss/len(omics2_train_loader):.4f}, RNA KL Loss: {total_rna_kl_loss/len(omics1_train_loader):.4f}, ATAC KL Loss: {total_atac_kl_loss/len(omics2_train_loader):.4f}, Acc_train: {ACC_Train:.2f}, Acc_test: {ACC_Test:.2f}')

    return E_rna, E_atac, D_rna, D_atac, Classifier


def train_single_omics(E_rna, D_rna, Classifier, optimizer_E_rna, optimizer_D_rna, optimizer_fc, scheduler_E_rna, scheduler_D_rna, scheduler_fc, omics1_train_loader, omics1_test_id_loader, num_epochs, device):
    E_rna.train()
    D_rna.train()
    Classifier.train()

    Loss_SupCon = SupConLoss(temperature=0.07)  # 0.07
    for epoch in range(num_epochs):
        total_loss_cls = 0.0
        total_rna_recon_loss = 0.0
        total_rna_kl_loss = 0.0
        for _, x_rna, y_rna in omics1_train_loader:
            x_rna = x_rna.to(device)

            y = y_rna.to(device)

            z_rna, mu_rna, logvar_rna = E_rna(x_rna)

            # ################################################
            # # 分类时，把双曲表示映射回切空间
            # z_h_tangent = logmap0(z_rna_h, c=1.0)
            #
            # # 融合欧氏特征和双曲切空间特征
            # z_rna_fusion = torch.cat([z_rna_e, z_h_tangent], dim=-1)
            # ################################################

            re_z_rna = reparameterize(mu_rna, logvar_rna)

            logits = Classifier(z_rna)

            loss_cls = F.cross_entropy(logits, y)

            # ---------------- decoder -----------------
            x_hat_rna = D_rna(re_z_rna)

            rna_recon_loss = F.mse_loss(x_hat_rna, x_rna, reduction="mean")

            rna_kl_loss = -0.5 * torch.mean(1 + logvar_rna - mu_rna.pow(2) - logvar_rna.exp())

            loss_sc_e = Loss_SupCon(F.normalize(z_rna), y) + Loss_SupCon(F.normalize(re_z_rna), y)

            loss = 0.1*loss_cls + rna_kl_loss + rna_recon_loss + 0.5*loss_sc_e  #

            # 清空梯度
            optimizer_E_rna.zero_grad()
            optimizer_D_rna.zero_grad()
            optimizer_fc.zero_grad()

            loss.backward()  # 反向传播

            # 更新权重
            optimizer_E_rna.step()
            optimizer_D_rna.step()
            optimizer_fc.step()

            # 更新学习率
            scheduler_E_rna.step()
            scheduler_D_rna.step()
            scheduler_fc.step()

            total_loss_cls += loss_cls
            total_rna_recon_loss += rna_recon_loss
            total_rna_kl_loss += rna_kl_loss

        ACC_Train = evaluate_model_with_single_omics(E_rna, Classifier, omics1_train_loader, device)

        if omics1_test_id_loader is not None:
            ACC_Test = evaluate_model_with_single_omics(E_rna, Classifier, omics1_test_id_loader, device)
        else:
            ACC_Test = ACC_Train

        # print(f'Epoch [{epoch+1}/{num_epochs}], CLS Loss: {total_loss_cls/len(omics1_train_loader):.4f}, RNA Recon Loss: {total_rna_recon_loss/len(omics1_train_loader):.4f}, RNA KL Loss: {total_rna_kl_loss/len(omics1_train_loader):.4f}, Acc_train: {ACC_Train:.2f}, Acc_test: {ACC_Test:.2f}')

        # if ACC_Train > 90.0:
        #     break

    return E_rna, D_rna, Classifier


def evaluate_model_with_multiomics(E_rna, E_atac, Classifier, paired_loader, device):
    E_rna.eval()
    E_atac.eval()
    Classifier.eval()

    correct = 0
    total = 0
    with torch.no_grad():  # 关闭梯度计算
        for _, x_rna, _, x_atac, labels in paired_loader:
            x_rna = x_rna.to(device)
            x_atac = x_atac.to(device)
            y = labels.to(device)

            z_rna, _, _ = E_rna(x_rna)
            z_atac, _, _ = E_atac(x_atac)

            z = torch.cat((z_rna, z_atac), dim=0)

            logits_all = Classifier(z)  # [2B, C]

            half = logits_all.size(0) // 2
            logits_r = logits_all[:half]
            logits_a = logits_all[half:]

            logits = (7*logits_r + 3*logits_a) / 10

            _, predicted = torch.max(logits.data, 1)

            total += logits.size(0)
            correct += (predicted == y).sum().item()
    accuracy = 100 * correct / total
    # print(f'Accuracy of the model on the test images: {accuracy:.2f}%')
    return accuracy


def get_pseudo_labeled_mixed_loader_with_single_omics(E_rna, Classifier, omics1_mixed_loader, threshold, device):
    E_rna.eval()
    Classifier.eval()

    omics1_data = []
    omics1_raw_data = []
    omics1_label = []

    with torch.no_grad():
        for x_raw_1, x_1, y_1 in omics1_mixed_loader:
            x_raw_1 = x_raw_1.to(device)
            x_1 = x_1.to(device)
            y_1 = y_1.to(device)
            z_1, _, _ = E_rna(x_1)

            output = Classifier(z_1)

            probabilities = F.softmax(output, dim=1)
            probabilities, predicted = torch.max(probabilities.data, 1)

            higher_confidence_mast = probabilities >= threshold
            omics1_data.append(x_1[higher_confidence_mast].detach().cpu())
            omics1_raw_data.append(x_raw_1[higher_confidence_mast].detach().cpu())
            omics1_label.append(predicted[higher_confidence_mast].detach().cpu())

            E_rna.zero_grad()
            Classifier.zero_grad()

    X = torch.cat(omics1_data, dim=0)
    label = torch.cat(omics1_label, dim=0)
    X_raw = torch.cat(omics1_raw_data, dim=0)

    filtered_omics1_dataset = BuildDataset(X_raw, X, label)

    return filtered_omics1_dataset


@torch.no_grad()
def get_latent_with_single_omics(
        E_rna,
        dataloader_rna,
        device,
        save_dir,
        omics=1,
        use_mu=True,  # True: mu, False: sampled z
        max_samples=2000,
        perplexity=30,
        random_state=0
):
    """
    获取低维表征Z，标签Y，以及数据X
    Joint t-SNE in latent z space for RNA / ATAC encoders
    """

    import os
    import numpy as np
    import torch
    import matplotlib.pyplot as plt
    from sklearn.manifold import TSNE

    os.makedirs(save_dir, exist_ok=True)

    E_rna.eval()

    # --------------------------------------------------
    # helper
    # --------------------------------------------------
    def collect_latent(E, dataloader, tag):
        z_list = []
        y_list = []
        x_list = []
        raw_x_list = []
        for batch in dataloader:
            if isinstance(batch, (list, tuple)):
                raw, x, y = batch
            else:
                x = batch

            x_list.append(x)
            raw_x_list.append(raw)

            x = x.to(device)
            z, mu, logvar = E(x)

            z_list.append(z.cpu())
            y_list.append(y.cpu())

            # if sum(z.shape[0] for z in z_list) >= max_samples:
            #     break

        Z = torch.cat(z_list, dim=0)  # [:max_samples]
        Y = torch.cat(y_list, dim=0)
        X = torch.cat(x_list, dim=0)
        raw_X = torch.cat(raw_x_list, dim=0)
        print(f"[{tag}] collected {Z.shape[0]} samples")
        return Z, Y, X, raw_X

    # --------------------------------------------------
    # collect latent representations
    # --------------------------------------------------
    Z_rna, Y_rna, X_rna, raw_X = collect_latent(E_rna, dataloader_rna, "RNA")

    return Z_rna, Y_rna, X_rna, raw_X


@torch.no_grad()
def visualize_decoder_tsne(
    E_rna,
    E_atac,
    D_rna,
    D_atac,
    dataloader,
    device,
    save_dir,
    modality="RNA",          # "RNA" or "ATAC"
    use_mu=True,             # True: 用 mu；False: 用 z 采样
    max_samples=2000,        # t-SNE 样本上限
    perplexity=30,
    random_state=0
):
    """
    Encode -> Decode -> Save -> t-SNE visualization

    Parameters
    ----------
    E_rna, E_atac : Encoder
    D_rna, D_atac : Decoder
    dataloader : DataLoader
        yield (x, y) or x
    device : torch.device
    save_dir : str
    modality : str
        "RNA" or "ATAC"
    use_mu : bool
        whether to use mu or reparameterized z
    """

    os.makedirs(save_dir, exist_ok=True)

    E = E_rna if modality == "RNA" else E_atac
    D = D_rna if modality == "RNA" else D_atac

    E.eval()
    D.eval()

    recon_list = []
    label_list = []

    for batch in dataloader:
        if isinstance(batch, (list, tuple)):
            _, x, y = batch
            label_list.append(y)
        else:
            x = batch

        x = x.to(device)

        embedding, mu, logvar = E(x)

        if use_mu:
            z = mu
        else:
            std = torch.exp(0.5 * logvar)
            z = mu + std * torch.randn_like(std)

        x_hat = D(z)

        recon_list.append(x_hat.cpu())

        if sum(r.shape[0] for r in recon_list) >= max_samples:
            break

    X_recon = torch.cat(recon_list, dim=0)[:max_samples]
    X_recon_np = X_recon.numpy()

    # 保存重建结果
    torch.save(X_recon, os.path.join(save_dir, f"{modality}_reconstruction.pt"))
    np.save(os.path.join(save_dir, f"{modality}_reconstruction.npy"), X_recon_np)

    if label_list:
        y_all = torch.cat(label_list, dim=0)[:max_samples].cpu().numpy()
    else:
        y_all = None

    # ---------- t-SNE ----------
    tsne = TSNE(
        n_components=2,
        perplexity=perplexity,
        init="pca",
        random_state=random_state
    )
    Z = tsne.fit_transform(X_recon_np)

    # ---------- plot ----------
    plt.figure(figsize=(6, 6))
    if y_all is not None:
        scatter = plt.scatter(Z[:, 0], Z[:, 1], c=y_all, s=6, cmap="tab20")
        plt.legend(*scatter.legend_elements(), title="Class", fontsize=8)
    else:
        plt.scatter(Z[:, 0], Z[:, 1], s=6)

    title_suffix = "mu" if use_mu else "z"
    plt.title(f"{modality} Decoder Reconstruction (t-SNE, {title_suffix})")
    plt.tight_layout()

    fig_path = os.path.join(
        save_dir, f"{modality}_decoder_tsne_{title_suffix}.png"
    )
    plt.savefig(fig_path, dpi=300)
    plt.close()

    print(f"[Saved]")
    print(f"  Recon tensor : {save_dir}/{modality}_reconstruction.pt")
    print(f"  t-SNE figure : {fig_path}")


def leiden_on_reconstruction_fixed(
    recon_x,
    save_dir,
    n_neighbors=15,
    n_pcs=50,
    resolution=1.0,
    var_threshold=1e-8,
    random_state=0
):
    import scanpy as sc
    import numpy as np
    import os
    import torch

    os.makedirs(save_dir, exist_ok=True)

    if isinstance(recon_x, torch.Tensor):
        X = recon_x.cpu().numpy()
    else:
        X = recon_x

    adata = sc.AnnData(X)

    # -------- 关键修复：过滤 near-zero variance features --------
    var = np.var(adata.X, axis=0)
    keep = var > var_threshold
    adata = adata[:, keep].copy()

    print(f"[Info] kept {keep.sum()} / {len(keep)} features")

    # -------- 正常流程 --------
    sc.pp.scale(adata)
    sc.tl.pca(adata, n_comps=n_pcs, svd_solver="arpack")

    sc.pp.neighbors(
        adata,
        n_neighbors=n_neighbors,
        n_pcs=n_pcs,
        random_state=random_state
    )
    sc.tl.leiden(adata, resolution=resolution)

    sc.tl.umap(adata, random_state=random_state)
    sc.pl.umap(adata, color="leiden", show=False, save="_leiden.png")

    adata.write(os.path.join(save_dir, "recon_leiden.h5ad"))
    return adata


import os
import torch
import numpy as np


@torch.no_grad()
def save_decoder_reconstruction(
    E_rna,
    E_atac,
    D_rna,
    D_atac,
    dataloader,
    device,
    save_dir,
    modality="RNA",          # "RNA" or "ATAC"
    use_mu=True,             # True: 用 mu；False: 用 z 采样
    max_samples=2000
):
    """
    Encode -> Decode -> Save reconstruction + labels

    Saves:
    ------
    {modality}_reconstruction.pt
    {modality}_reconstruction.npy
    {modality}_labels.pt (if available)
    {modality}_labels.npy (if available)
    """

    os.makedirs(save_dir, exist_ok=True)

    E = E_rna if modality == "RNA" else E_atac
    # D = D_rna if modality == "RNA" else D_atac  # 重建
    D = D_atac if modality == "RNA" else D_rna  # 预测

    E.eval()
    D.eval()

    recon_list = []
    label_list = []

    n_collected = 0

    for batch in dataloader:
        if isinstance(batch, (list, tuple)):
            # 按你现在的 dataloader 结构
            _, x, y = batch
            label_list.append(y)
        else:
            x = batch
            y = None

        x = x.to(device)

        embedding, mu, logvar = E(x)

        if use_mu:
            z = mu
        else:
            std = torch.exp(0.5 * logvar)
            z = mu + std * torch.randn_like(std)

        x_hat = D(z)

        recon_list.append(x_hat.cpu())
        n_collected += x_hat.shape[0]

        if n_collected >= max_samples:
            break

    # ---------- concat & trim ----------
    X_recon = torch.cat(recon_list, dim=0)[:max_samples]

    # ---------- save reconstruction ----------
    torch.save(
        X_recon,
        os.path.join(save_dir, f"{modality}_reconstruction.pt")
    )
    np.save(
        os.path.join(save_dir, f"{modality}_reconstruction.npy"),
        X_recon.numpy()
    )

    # ---------- save labels (if any) ----------
    if label_list:
        y_all = torch.cat(label_list, dim=0)[:max_samples]
        torch.save(
            y_all,
            os.path.join(save_dir, f"{modality}_labels.pt")
        )
        np.save(
            os.path.join(save_dir, f"{modality}_labels.npy"),
            y_all.cpu().numpy()
        )

    print(f"[Saved reconstruction]")
    print(f"  X shape : {X_recon.shape}")
    print(f"  Path    : {save_dir}/{modality}_reconstruction.pt")

    if label_list:
        print(f"[Saved labels]")
        print(f"  y shape : {y_all.shape}")
        print(f"  Path    : {save_dir}/{modality}_labels.pt")

    return X_recon, (y_all if label_list else None)












############################################
###############     ADBC     ###############
############################################


















############################################
##############     DNCTE     ###############
############################################