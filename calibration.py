import torch
import numpy as np
import torch.nn.functional as F
from display_results import get_measures, print_measures
from data_split import BuildDataset
from torch.utils.data import DataLoader
import random
from collections import defaultdict
from algorithm_utils import evaluate_model
import scanpy as sc
from sklearn.metrics.pairwise import cosine_similarity
from scDiscovery_architecture import logmap0
from threshold_estimate import get_scores_and_density_plot, threshold_msp_simple, plot_roc_with_threshold


seed = 0
random.seed(0)
torch.manual_seed(seed)
np.random.seed(seed)
torch.cuda.manual_seed(seed)
rng = np.random.default_rng(seed)
torch.cuda.manual_seed_all(seed)

concat = lambda x: np.concatenate(x, axis=0)
to_np = lambda x: x.data.cpu().numpy()


def pseudo_labelling_for_test_mixed_data(X_mixed_test, Y_mixed_test, Z_mixed_test, X_mixed_test_raw, n_neighbors=15, resolution=0.45, omics='RNA'):
    """
    在低维表征空间Z中聚类，打伪标签
    """

    adata = sc.AnnData(Z_mixed_test.numpy())

    sc.pp.neighbors(adata, n_neighbors=n_neighbors, metric="euclidean")
    sc.tl.leiden(adata, resolution=resolution)

    cluster_labels = adata.obs["leiden"]

    adata.obsm["X_mixed_test"] = X_mixed_test.numpy()
    adata.obsm["X_mixed_test_raw"] = X_mixed_test_raw.numpy()

    adata.obs['label'] = Y_mixed_test
    adata.obs["label"] = adata.obs["label"].astype("category")

    sc.tl.umap(adata)
    sc.tl.tsne(adata)
    sc.pl.umap(adata, color=["leiden", "label"],
               frameon=False,
               legend_loc="on data",
               ncols=2,
               wspace=0.5,
               save=f"{omics}_leiden.png"
               )
    return adata


# Step 1：计算已知 8 类的质心
def compute_class_centroids(Z, Y, num_classes=8):
    centroids = []
    for c in range(num_classes):
        centroids.append(Z[Y == c].mean(axis=0))
    return np.stack(centroids)  # [8, D]


# Step 2：计算混合数据每个 Leiden 簇的质心
def compute_cluster_centroids(Z_mixed, cluster_labels):
    clusters = np.unique(cluster_labels)
    centroids = {}
    for c in clusters:
        centroids[c] = Z_mixed[cluster_labels == c].mean(axis=0)
    return centroids  # dict: {cluster_id: centroid}


def distinguish_known_unknown(Z_train, Y_train, mixed_test_adata, known_class_number):
    """
    区分已知类和未知类
    """

    # Step 1：计算已知 8 类的质心
    known_centroids = compute_class_centroids(Z_train, Y_train, num_classes=known_class_number)

    # 取特征
    Z_mixed = mixed_test_adata.obsm["Z"] if "Z" in mixed_test_adata.obsm else mixed_test_adata.X
    cluster_labels = mixed_test_adata.obs["leiden"].astype(int).values

    # Step 2：计算混合数据每个 Leiden 簇的质心
    cluster_centroids = compute_cluster_centroids(Z_mixed, cluster_labels)

    # Step 3：计算“簇 ↔ 已知类”的相似度对齐
    cluster_ids = list(cluster_centroids.keys())
    C = np.stack([cluster_centroids[c] for c in cluster_ids])  # [K, D]

    # [K, 8]
    sim = cosine_similarity(C, known_centroids)

    # 每个簇最像哪个已知类
    best_class = sim.argmax(axis=1)
    best_score = sim.max(axis=1)

    # Step 4：为 8 个已知类选最匹配的 8 个簇
    selected_clusters = set()

    for k in range(known_class_number):
        idx = np.argmax(sim[:, k])  # 与第 k 类最像的簇
        selected_clusters.add(cluster_ids[idx])

    known_like_clusters = sorted(list(selected_clusters))

    # Step 5：其余簇 = 未知类簇
    all_clusters = set(cluster_ids)
    unknown_clusters = sorted(list(all_clusters - set(known_like_clusters)))

    # Step 6：给每个细胞打“已知 / 未知”标签
    def assign_known_unknown(cluster_labels, known_clusters):
        flags = []
        for c in cluster_labels:
            flags.append("known_like" if str(c) in known_clusters else "unknown")
        return flags

    mixed_test_adata.obs["known_unknown"] = assign_known_unknown(
        mixed_test_adata.obs["leiden"].astype(str).values,
        set(map(str, known_like_clusters))
    )

    sc.pl.umap(mixed_test_adata, color=["leiden", "known_unknown"], frameon=False, legend_loc='on data', save='known_unknown.png')

    return mixed_test_adata


def get_and_remapping_novel_cell_type(mixed_test_adata, known_class_number=8):
    """
    获取未知类细胞，并映射伪标签，接在已知类的标签后边
    """

    # Step 1：取出 unknown 细胞
    unknown_mask = mixed_test_adata.obs["known_unknown"] == "unknown"
    adata_unknown = mixed_test_adata[unknown_mask].copy()

    # Step 2：获取 unknown 簇列表
    unknown_clusters = sorted(
        adata_unknown.obs["leiden"].astype(str).unique()
    )

    # Step 3：建立“旧簇 → 新标签”的映射
    start_label = known_class_number  # 已知类是 0~7

    cluster2new = {
        c: start_label + i
        for i, c in enumerate(unknown_clusters)
    }

    # Step 4：给 unknown 细胞打新标签
    adata_unknown.obs["new_label"] = (
        adata_unknown.obs["leiden"].astype(str).map(cluster2new).astype(int)
    )

    return adata_unknown


def get_ood_scores(E_rna, E_atac, Classifier, omics1_loader, omics2_loader, device):
    E_rna.eval()
    E_atac.eval()
    Classifier.eval()

    omics1_score = []
    omics2_score = []

    omics1_data = []
    omics2_data = []

    omics1_raw_data = []
    omics2_raw_data = []

    omics1_label = []
    omics2_label = []

    with torch.no_grad():
        for (x_raw_1, x_1, y_1), (x_raw_2, x_2, y_2) in zip(omics1_loader, omics2_loader):
            z_1, _, _ = E_rna(x_1.to(device))
            z_2, _, _ = E_atac(x_2.to(device))
            z = torch.cat((z_1, z_2), 0)

            output = Classifier(z)

            smax = to_np(F.softmax(output, dim=1))

            scores = np.max(smax, axis=1)  # np.max(smax, axis=1)

            omics1_data.append(x_1)
            omics2_data.append(x_2)

            omics1_raw_data.append(x_raw_1)
            omics2_raw_data.append(x_raw_2)

            omics1_label.append(y_1)
            omics2_label.append(y_2)

            omics1_score.append(scores[:len(x_1)])
            omics2_score.append(scores[len(x_1):])

            E_rna.zero_grad()
            E_atac.zero_grad()
            Classifier.zero_grad()

    return concat(omics1_score).copy(), concat(omics2_score).copy(), torch.cat(omics1_data, dim=0), torch.cat(omics2_data, dim=0), torch.cat(omics1_label, dim=0), torch.cat(omics2_label, dim=0), torch.cat(omics1_raw_data, dim=0), torch.cat(omics2_raw_data, dim=0)


def get_ood_scores_with_single_omics(E_rna, Classifier, omics1_loader, device):
    E_rna.eval()
    Classifier.eval()

    omics1_score = []
    omics1_data = []
    omics1_raw_data = []
    omics1_label = []

    with torch.no_grad():
        for x_raw_1, x_1, y_1 in omics1_loader:
            z_1, _, _ = E_rna(x_1.to(device))

            output = Classifier(z_1)

            smax = to_np(F.softmax(output, dim=1))

            scores = np.max(smax, axis=1)  # np.max(smax, axis=1)

            omics1_data.append(x_1)

            omics1_raw_data.append(x_raw_1)

            omics1_label.append(y_1)

            omics1_score.append(scores[:len(x_1)])

            E_rna.zero_grad()
            Classifier.zero_grad()

    return concat(omics1_score).copy(), torch.cat(omics1_data, dim=0), torch.cat(omics1_label, dim=0), torch.cat(omics1_raw_data, dim=0)


def get_ood_scores_with_multiomics(E_rna, E_atac, Classifier, paired_loader, device):
    E_rna.eval()
    Classifier.eval()

    paired_score = []
    omics_rna_data = []
    omics_rna_raw_data = []
    omics_atac_score = []
    omics_atac_data = []
    omics_atac_raw_data = []
    paired_label = []

    with torch.no_grad():
        for x_raw_rna, x_rna, x_raw_atac, x_atac, labels in paired_loader:
            z_rna, _, _ = E_rna(x_rna.to(device))
            z_atac, _, _ = E_atac(x_atac.to(device))

            z = torch.cat((z_rna, z_atac), dim=0)

            output_all = Classifier(z)  # [2B, C]

            half = output_all.size(0) // 2
            logits_r = output_all[:half]
            logits_a = output_all[half:]

            logits = (7*logits_r + 3*logits_a) / 10

            smax = to_np(F.softmax(logits, dim=1))

            scores = np.max(smax, axis=1)  # np.max(smax, axis=1)

            omics_rna_data.append(x_rna)
            omics_atac_data.append(x_atac)

            omics_rna_raw_data.append(x_raw_rna)
            omics_atac_raw_data.append(x_raw_atac)

            paired_label.append(labels)

            paired_score.append(scores[:len(x_rna)])

            E_rna.zero_grad()
            E_atac.zero_grad()
            Classifier.zero_grad()

    return concat(paired_score).copy(), torch.cat(omics_rna_data, dim=0), torch.cat(omics_atac_data, dim=0), torch.cat(paired_label, dim=0), torch.cat(omics_rna_raw_data, dim=0), torch.cat(omics_atac_raw_data, dim=0)


@torch.no_grad()
def compute_class_centers(
    encoder,
    train_loader,
    device="cuda",
):
    """
    encoder: E_rna 或 E_atac
    train_loader: omics1_train_loader 或 omics2_train_loader
                  batch = (x, y)
    return:
      centers: dict {class_id: center_tensor [D]}
    """
    encoder.eval()

    feats_by_class = defaultdict(list)

    for _, x, y in train_loader:
        x = x.to(device)
        y = y.to(device)

        z, _, _ = encoder(x)   # [B, D]

        for i in range(len(y)):
            c = y[i].item()
            feats_by_class[c].append(z[i])

    centers = {
        c: torch.stack(feats).mean(dim=0)
        for c, feats in feats_by_class.items()
    }

    return centers


@torch.no_grad()
def select_farthest_data(
    encoder,
    mixed_loader,
    centers_dict,
    ratio=0.5,
    device="cuda"
):
    """
    不依赖 dataset index，不怕 shuffle=True
    返回：selected_x, selected_y （可选 selected_dist）
    """
    encoder.eval()
    centers = torch.stack(list(centers_dict.values())).to(device)  # [C, D]

    xs, ys, ds = [], [], []

    for _, x, y in mixed_loader:
        x = x.to(device)
        z, _, _ = encoder(x)                      # [B, D]
        dist = torch.cdist(z, centers)      # [B, C]
        min_dist = dist.min(dim=1)[0]       # [B]

        xs.append(x.detach().cpu())
        ys.append(y.detach().cpu())
        ds.append(min_dist.detach().cpu())

    X = torch.cat(xs, dim=0)   # [N, ...]
    Y = torch.cat(ys, dim=0)   # [N]
    D = torch.cat(ds, dim=0)   # [N]

    k = int(len(D) * ratio)
    sel = torch.argsort(D, descending=True)[:k]

    data = X[sel]
    labels = Y[sel]

    filtered_dataset = BuildDataset(labels, data, labels)

    filtered_loader = DataLoader(
        dataset=filtered_dataset, batch_size=128, shuffle=True, drop_last=False
    )

    return filtered_loader


def init_filter(E_rna, E_atac, omics1_train_loader, omics2_train_loader, omics1_mixed_loader, omics2_mixed_loader, omics1_test_ood_loader, omics2_test_ood_loader):
    # 1) 计算类质心（train 数据）
    omics1_center = compute_class_centers(
        E_rna, omics1_train_loader
    )

    omics2_center = compute_class_centers(
        E_atac, omics2_train_loader
    )

    # 2) mixed 数据 → 距离
    init_omics1_filtered_loader = select_farthest_data(
        E_rna,
        omics1_mixed_loader,
        omics1_center,
        ratio=0.4,
        device="cuda"
    )

    init_omics2_filtered_loader = select_farthest_data(
        E_atac,
        omics2_mixed_loader,
        omics2_center,
        ratio=0.4,
        device="cuda"
    )

    return init_omics1_filtered_loader, init_omics2_filtered_loader


def init_filter_with_single_omics(E_rna, omics1_train_loader, omics1_mixed_loader, omics1_test_ood_loader):
    # 1) 计算类质心（train 数据）
    omics1_center = compute_class_centers(
        E_rna, omics1_train_loader
    )

    # 2) mixed 数据 → 距离
    init_omics1_filtered_loader = select_farthest_data(
        E_rna,
        omics1_mixed_loader,
        omics1_center,
        ratio=0.2,
        device="cuda"
    )

    return init_omics1_filtered_loader


def filter_ood_data(E_rna, E_atac, Classifier, omics1_train_loader, omics2_train_loader, omics1_mixed_loader, omics2_mixed_loader, omics1_test_ood_loader, omics2_test_ood_loader, filtered_omics1_ood_loader, filtered_omics2_ood_loader,class_number, pre_threshold, num_iters, iter, recall_level=0.95):
    id_omics1_scores, id_omics2_scores, _, _, _, _, _, _ = get_ood_scores(E_rna, E_atac, Classifier, omics1_train_loader, omics2_train_loader)
    id_scores = np.concatenate((id_omics1_scores, id_omics2_scores), axis=0)

    mixed_omics1_scores, mixed_omics2_scores, mixed_omics1_data, mixed_omics2_data, _, _, mixed_omics1_raw_data, mixed_omics2_raw_data = get_ood_scores(E_rna, E_atac,
                                                                                                    Classifier,
                                                                                                    omics1_mixed_loader,
                                                                                                    omics2_mixed_loader)  # omics1_mixed_loader, omics2_mixed_loader


    measures = get_measures(id_scores, id_scores, recall_level=recall_level, plot=False)

    threshold = measures[3]
    if threshold > 0.3:
        threshold = 2*(1/class_number)

    print(f'Threshold in This Iter:{threshold}')

    mask1 = mixed_omics1_scores < threshold
    mask2 = mixed_omics2_scores < threshold
    filtered_data1 = mixed_omics1_data[mask1]  # 筛选后的数据
    filtered_data2 = mixed_omics2_data[mask2]  # 筛选后的数据

    filtered_raw_data1 = mixed_omics1_raw_data[mask1]
    filtered_raw_data2 = mixed_omics2_raw_data[mask2]

    if len(filtered_data1) == 0 or len(filtered_data2) == 0:
        threshold += 0.1

        mask1 = mixed_omics1_scores < threshold
        mask2 = mixed_omics2_scores < threshold
        filtered_data1 = mixed_omics1_data[mask1]  # 筛选后的数据
        filtered_data2 = mixed_omics2_data[mask2]  # 筛选后的数据

        filtered_raw_data1 = mixed_omics1_raw_data[mask1]
        filtered_raw_data2 = mixed_omics2_raw_data[mask2]

    filtered_label1 = 1000*torch.ones(len(filtered_data1), dtype=torch.long)
    filtered_label2 = 1000*torch.ones(len(filtered_data2), dtype=torch.long)

    filtered_omics1_dataset = BuildDataset(filtered_raw_data1, filtered_data1, filtered_label1)
    filtered_omics2_dataset = BuildDataset(filtered_raw_data2, filtered_data2, filtered_label2)

    filtered_omics1_loader = DataLoader(
        dataset=filtered_omics1_dataset, batch_size=128, shuffle=True, drop_last=False
    )

    filtered_omics2_loader = DataLoader(
        dataset=filtered_omics2_dataset, batch_size=128, shuffle=True, drop_last=False
    )

    return filtered_omics1_loader, filtered_omics2_loader, threshold


def filter_ood_data_with_single_omics(E_rna, Classifier, omics1_train_loader, omics1_mixed_loader, omics1_test_ood_loader, filtered_omics1_ood_loader, class_number, pre_threshold, num_iters, device, iter, recall_level=0.95):
    id_scores, _, _, _ = get_ood_scores_with_single_omics(E_rna, Classifier, omics1_train_loader, device)

    mixed_omics1_scores, mixed_omics1_data, _, mixed_omics1_raw_data = get_ood_scores_with_single_omics(E_rna, Classifier,
                                                                                                    omics1_mixed_loader, device)  # omics1_mixed_loader, omics2_mixed_loader

    # if filtered_omics1_ood_loader == None:
    #     mixed_scores = np.concatenate((mixed_omics1_scores, mixed_omics2_scores), axis=0)
    # else:
    #     filtered_omics1_scores, filtered_omics2_scores, _, _, _, _ = get_ood_scores(E_rna, E_atac, Classifier, filtered_omics1_ood_loader, filtered_omics2_ood_loader)
    #     mixed_scores = np.concatenate((filtered_omics1_scores, filtered_omics2_scores), axis=0)

    # measures = get_measures(mixed_scores, id_scores, recall_level=recall_level, plot=False)
    measures = get_measures(id_scores, id_scores, recall_level=recall_level, plot=False)
    # print_measures(measures[0], measures[1], measures[2])

    threshold = pre_threshold - 0.0
    # threshold = measures[3]
    # if threshold > 0.3:
    #     threshold = 2*(1/class_number)

    # if pre_threshold is not None:
    #     threshold = pre_threshold-0.08
    # if iter == 1:
    #     threshold = 0.3
    # elif iter == 2:
    #     threshold = 0.15
    # elif iter == 3:
    #     threshold = 0.15
    # elif iter == 4:
    #     threshold = 0.1
    # elif iter == 5:
    #     threshold = 0.1

    # print(f'Threshold in This Iter:{threshold}')

    mask1 = mixed_omics1_scores < threshold
    filtered_data1 = mixed_omics1_data[mask1]  # 筛选后的数据

    filtered_raw_data1 = mixed_omics1_raw_data[mask1]

    if len(filtered_data1) == 0:
        # threshold += 0.1

        mask1 = mixed_omics1_scores < threshold+0.1
        filtered_data1 = mixed_omics1_data[mask1]  # 筛选后的数据

        filtered_raw_data1 = mixed_omics1_raw_data[mask1]

    filtered_label1 = 1000*torch.ones(len(filtered_data1), dtype=torch.long)

    filtered_omics1_dataset = BuildDataset(filtered_raw_data1, filtered_data1, filtered_label1)

    filtered_omics1_loader = DataLoader(
        dataset=filtered_omics1_dataset, batch_size=128, shuffle=True, drop_last=False
    )

    # print(f'Threshold in This Iter:{threshold}')

    return filtered_omics1_loader, threshold


def train_distinguish_id_ood_with_single_omics(E_rna, Classifier, optimizer_E_rna, optimizer_fc, scheduler_E_rna, scheduler_fc, omics1_train_loader,
        omics1_test_id_loader, omics1_test_ood_loader, omics1_mixed_loader, class_number, num_iters, num_epochs, recall_level, threshold, device):
    filtered_omics1_ood_loader = None
    # threshold = None

    if omics1_test_id_loader is None:
        omics1_test_id_loader = omics1_train_loader

    final_threshold = evaluate_ood_detection_with_single_omics(E_rna, Classifier, omics1_test_id_loader, omics1_test_ood_loader, recall_level, device)

    threshold = threshold
    for i in range(num_iters):
        # print(f'Iteration {i + 1}')
        if i == 0:
            filtered_omics1_ood_loader, threshold = filter_ood_data_with_single_omics(E_rna, Classifier,
                                                                                      omics1_train_loader,
                                                                                      omics1_mixed_loader,
                                                                                      omics1_test_ood_loader,
                                                                                      filtered_omics1_ood_loader,
                                                                                      class_number,
                                                                                      threshold,
                                                                                      num_iters,
                                                                                      device,
                                                                                      iter=i,
                                                                                      recall_level=0.95,
                                                                                      )
        else:
            filtered_omics1_ood_loader, threshold = filter_ood_data_with_single_omics(E_rna, Classifier, omics1_train_loader,
                                                                                     omics1_mixed_loader,
                                                                                     omics1_test_ood_loader,
                                                                                     filtered_omics1_ood_loader,
                                                                                     class_number,
                                                                                     threshold,
                                                                                     num_iters,
                                                                                     device,
                                                                                     iter=i,
                                                                                     recall_level=0.95)
        if i == 0:
            epochs = num_epochs
        else:
            epochs = num_epochs
        for epoch in range(epochs):
            E_rna.train()
            Classifier.train()
            for (_, x_id_1, y_1), (_, x_ood_1, _) in zip(omics1_train_loader, filtered_omics1_ood_loader):  # filtered_omics1_ood_loader omics1_test_ood_loader
                x_id_1, x_ood_1 = x_id_1.to(device), x_ood_1.to(device)
                y_1 = y_1.to(device)
                y = y_1

                x_1 = torch.cat((x_id_1, x_ood_1), 0)

                z_1, _, _ = E_rna(x_1)

                logits = Classifier(z_1)

                logits_1 = logits[:len(x_1)]

                logits_id_1 = logits_1[:len(x_id_1)]

                logits_ood_1 = logits_1[len(x_id_1):]

                l_ce = F.cross_entropy(logits_id_1, y)

                l_oe = 0.5 * -(logits_ood_1.mean(1) - torch.logsumexp(logits_ood_1, dim=1)).mean()

                loss = l_ce + l_oe

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

            final_threshold = evaluate_ood_detection_with_single_omics(E_rna, Classifier, omics1_test_id_loader, omics1_test_ood_loader, recall_level=0.95, device=device)

        # threshold -= 0.2
    threshold, filtered_omics1_ood_loader, filtered_omics1_id_loader = get_detected_ood_data_with_single_omics(E_rna, Classifier, omics1_mixed_loader, threshold, device=device)  # final_threshold

    return E_rna, Classifier, filtered_omics1_ood_loader, filtered_omics1_id_loader


def get_detected_ood_data(E_rna, E_atac, Classifier, omics1_mixed_loader, omics2_mixed_loader, final_threshold, device):
    E_rna.eval()
    E_atac.eval()
    Classifier.eval()

    mixed_omics1_scores, mixed_omics2_scores, mixed_omics1_data, mixed_omics2_data, mixed_omics1_label, mixed_omics2_label, mixed_omics1_raw_data, mixed_omics2_raw_data = get_ood_scores(E_rna, E_atac,
                                                                                                    Classifier,
                                                                                                    omics1_mixed_loader,
                                                                                                    omics2_mixed_loader, device)

    mask1_ood = mixed_omics1_scores < final_threshold
    mask2_ood = mixed_omics2_scores < final_threshold

    filtered_ood_data1 = mixed_omics1_data[mask1_ood]  # 筛选后的数据
    filtered_ood_data2 = mixed_omics2_data[mask2_ood]  # 筛选后的数据
    filtered_ood_label1 = mixed_omics1_label[mask1_ood]  # 筛选后的数据
    filtered_ood_label2 = mixed_omics2_label[mask2_ood]  # 筛选后的数据
    filtered_ood_raw_data1 = mixed_omics1_raw_data[mask1_ood]
    filtered_ood_raw_data2 = mixed_omics2_raw_data[mask2_ood]


    filtered_id_data1 = mixed_omics1_data[~mask1_ood]  # 筛选后的数据
    filtered_id_data2 = mixed_omics2_data[~mask2_ood]  # 筛选后的数据
    filtered_id_label1 = mixed_omics1_label[~mask1_ood]  # 筛选后的数据
    filtered_id_label2 = mixed_omics2_label[~mask2_ood]  # 筛选后的数据
    filtered_id_raw_data1 = mixed_omics1_raw_data[~mask1_ood]
    filtered_id_raw_data2 = mixed_omics2_raw_data[~mask2_ood]

    filtered_omics1_ood_dataset = BuildDataset(filtered_ood_raw_data1, filtered_ood_data1, filtered_ood_label1)
    filtered_omics2_ood_dataset = BuildDataset(filtered_ood_raw_data2, filtered_ood_data2, filtered_ood_label2)

    filtered_omics1_id_dataset = BuildDataset(filtered_id_raw_data1, filtered_id_data1, filtered_id_label1)
    filtered_omics2_id_dataset = BuildDataset(filtered_id_raw_data2, filtered_id_data2, filtered_id_label2)

    filtered_omics1_ood_loader = DataLoader(
        dataset=filtered_omics1_ood_dataset, batch_size=128, shuffle=True, drop_last=False
    )

    filtered_omics2_ood_loader = DataLoader(
        dataset=filtered_omics2_ood_dataset, batch_size=128, shuffle=True, drop_last=False
    )

    filtered_omics1_id_loader = DataLoader(
        dataset=filtered_omics1_id_dataset, batch_size=128, shuffle=True, drop_last=False
    )

    filtered_omics2_id_loader = DataLoader(
        dataset=filtered_omics2_id_dataset, batch_size=128, shuffle=True, drop_last=False
    )

    return filtered_omics1_ood_loader, filtered_omics2_ood_loader, filtered_omics1_id_loader, filtered_omics2_id_loader


def get_detected_ood_data_with_single_omics(
    E_rna,
    Classifier,
    omics1_mixed_loader,
    final_threshold,
    device,
    step=0.1,
    max_retry=10
):
    E_rna.eval()
    Classifier.eval()

    mixed_omics1_scores, mixed_omics1_data, mixed_omics1_label, mixed_omics1_raw_data = get_ood_scores_with_single_omics(
        E_rna,
        Classifier,
        omics1_mixed_loader,
        device
    )

    retry_count = 0

    while True:
        mask1_ood = mixed_omics1_scores < final_threshold

        filtered_ood_data1 = mixed_omics1_data[mask1_ood]
        filtered_ood_label1 = mixed_omics1_label[mask1_ood]
        filtered_ood_raw_data1 = mixed_omics1_raw_data[mask1_ood]

        filtered_id_data1 = mixed_omics1_data[~mask1_ood]
        filtered_id_label1 = mixed_omics1_label[~mask1_ood]
        filtered_id_raw_data1 = mixed_omics1_raw_data[~mask1_ood]

        filtered_omics1_ood_dataset = BuildDataset(
            filtered_ood_raw_data1,
            filtered_ood_data1,
            filtered_ood_label1
        )

        filtered_omics1_id_dataset = BuildDataset(
            filtered_id_raw_data1,
            filtered_id_data1,
            filtered_id_label1
        )

        try:
            filtered_omics1_ood_loader = DataLoader(
                dataset=filtered_omics1_ood_dataset,
                batch_size=128,
                shuffle=True,
                drop_last=False
            )

            filtered_omics1_id_loader = DataLoader(
                dataset=filtered_omics1_id_dataset,
                batch_size=128,
                shuffle=True,
                drop_last=False
            )

            # print(f"Final threshold used: {final_threshold}")
            # print(f"Detected OOD samples: {len(filtered_omics1_ood_dataset)}")
            # print(f"Detected ID samples: {len(filtered_omics1_id_dataset)}")

            return final_threshold, filtered_omics1_ood_loader, filtered_omics1_id_loader

        except ValueError as e:
            if "num_samples should be a positive integer" in str(e):
                retry_count += 1

                if retry_count > max_retry:
                    raise ValueError(
                        f"After {max_retry} retries, DataLoader still failed. "
                        f"Current threshold={final_threshold}, "
                        f"OOD samples={len(filtered_omics1_ood_dataset)}, "
                        f"ID samples={len(filtered_omics1_id_dataset)}"
                    )

                final_threshold += step
                print(
                    f"Empty dataset detected. "
                    f"Increase threshold to {final_threshold:.4f} "
                    f"and retry... [{retry_count}/{max_retry}]"
                )

            else:
                raise e


def evaluate_ood_detection(E_rna, E_atac, Classifier, omics1_test_id_loader, omics2_test_id_loader, omics1_test_ood_loader, omics2_test_ood_loader, recall_level, device):
    id_omics1_scores, id_omics2_scores, _, _, _, _, _, _ = get_ood_scores(E_rna, E_atac, Classifier, omics1_test_id_loader,
                                                              omics2_test_id_loader, device)
    ood_omics1_scores, ood_omics2_scores, _, _, _, _, _, _ = get_ood_scores(E_rna, E_atac, Classifier, omics1_test_ood_loader, omics2_test_ood_loader, device)

    id_scores = np.concatenate((id_omics1_scores, id_omics2_scores), axis=0)
    ood_scores = np.concatenate((ood_omics1_scores, ood_omics2_scores), axis=0)

    # measures = get_measures(ood_scores, id_scores, recall_level=recall_level, plot=False)
    measures = get_measures(id_scores, ood_scores, recall_level=recall_level, plot=False)
    print_measures(measures[0], measures[1], measures[2])

    print(measures[3])

    return measures[3]


def evaluate_ood_detection_with_single_omics(E_rna, Classifier, omics1_test_id_loader, omics1_test_ood_loader, recall_level, device):
    id_scores, _, _, _ = get_ood_scores_with_single_omics(E_rna, Classifier, omics1_test_id_loader, device)
    ood_scores, _, _, _ = get_ood_scores_with_single_omics(E_rna, Classifier, omics1_test_ood_loader, device)

    # measures = get_measures(ood_scores, id_scores, recall_level=recall_level, plot=False)
    measures = get_measures(id_scores, ood_scores, recall_level=recall_level, plot=False)
    # print_measures(measures[0], measures[1], measures[2])

    # print(measures[3])

    return measures[3]


def adaptive_decision_boundary_calibration(E_rna, Classifier, optimizer_E_rna, optimizer_fc, scheduler_E_rna, scheduler_fc, omics1_train_loader,
        omics1_test_id_loader, omics1_test_ood_loader, omics1_mixed_loader, class_number, device, dataset_name, setting, num_iters=1,
        num_epochs=5, recall_level=0.95):
    print("=================== Before Calibration ===================")
    scores, id_scores, ood_scores = get_scores_and_density_plot(E_rna, Classifier, omics1_test_id_loader, omics1_test_ood_loader,
                                         omics1_train_loader, dataset_name, setting, device)

    threshold_before, _, _ = threshold_msp_simple(scores, mode_offset=0.1)

    E_rna, Classifier, filtered_omics1_ood_loader, filtered_omics1_id_loader = train_distinguish_id_ood_with_single_omics(
        E_rna, Classifier, optimizer_E_rna, optimizer_fc, scheduler_E_rna, scheduler_fc, omics1_train_loader,
        omics1_test_id_loader, omics1_test_ood_loader, omics1_mixed_loader, class_number=class_number, num_iters=num_iters,
        num_epochs=num_epochs, recall_level=recall_level, threshold=threshold_before, device=device)

    print("=================== After Calibration ===================")

    scores, id_scores, ood_scores = get_scores_and_density_plot(E_rna, Classifier, omics1_test_id_loader, omics1_test_ood_loader,
                                         omics1_train_loader, dataset_name, setting, device)

    threshold_after, _, _ = threshold_msp_simple(scores, mode_offset=0.1)

    if threshold_after != threshold_before:
        threshold, filtered_omics1_ood_loader, filtered_omics1_id_loader = get_detected_ood_data_with_single_omics(E_rna, Classifier,
                                                                                                    omics1_mixed_loader,
                                                                                                    threshold_after,
                                                                                                    device=device)
    else:
        threshold = threshold_before

    print("=================== ROC Curve ===================")

    result = plot_roc_with_threshold(
        id_scores,
        ood_scores,
        threshold=threshold,
        save_path=f"./figures/roc_with_threshold_{dataset_name}_{setting}_trained.pdf"
    )

    # print(result)

    return threshold, filtered_omics1_ood_loader, filtered_omics1_id_loader