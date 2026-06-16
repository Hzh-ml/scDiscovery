import torch
import numpy as np
import scanpy as sc
import anndata as ad
import pandas as pd
import torch.nn.functional as F
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score, normalized_mutual_info_score, homogeneity_score
from torch.utils.data import TensorDataset, DataLoader
from scipy.sparse import csr_matrix
import random
import warnings
warnings.filterwarnings("ignore")
from sklearn.feature_extraction.text import TfidfTransformer
from data_split import BuildDataset, DataLoader
from contrastive import expand_classifier
import torch.optim as optim
import matplotlib.pyplot as plt
from algorithm_utils import train_single_omics

seed = 0
random.seed(0)
torch.manual_seed(seed)
np.random.seed(seed)
torch.cuda.manual_seed(seed)
rng = np.random.default_rng(seed)
torch.cuda.manual_seed_all(seed)

to_np = lambda x: x.data.cpu().numpy()


def pseudo_label(
        dataloader,
        n_pcs=30,
        n_neighbors=15,
        resolution=0.5,
        use_scale=False,
        device="cpu",
        data_type="rna",  # 'rna'  'atac'
        atac_n_top_features=None,
):
    """
    Generate pseudo labels for OOD data using dimensionality reduction + clustering.
    Stores true labels from dataloader in adata.obs for comparison.
    """
    # --------------------------------------------------
    # 1. Collect data AND true labels
    # --------------------------------------------------
    X_list = []
    true_labels_list = []

    for batch_data in dataloader:
        if isinstance(batch_data, (list, tuple)) and len(batch_data) >= 2:
            x, y = batch_data[0], batch_data[1]
        else:
            x = batch_data
            y = None

        X_list.append(x.detach().cpu())

        if y is not None:
            if isinstance(y, torch.Tensor):
                true_labels_list.append(y.detach().cpu().numpy())
            else:
                true_labels_list.append(y)

    X = torch.cat(X_list, dim=0).numpy()  # [N, D]

    adata = ad.AnnData(X)

    if true_labels_list:
        if all(isinstance(arr, np.ndarray) for arr in true_labels_list):
            true_labels = np.concatenate(true_labels_list, axis=0)
        else:
            true_labels = np.array(true_labels_list).flatten()

        if len(true_labels) == adata.n_obs:
            adata.obs["true_label"] = true_labels
            adata.obs["true_label"] = adata.obs["true_label"].astype('category')
            print(f"Added true labels to adata.obs['true_label'], shape: {true_labels.shape}")
        else:
            print(f"Warning: True labels count ({len(true_labels)}) doesn't match cell count ({adata.n_obs})")

    # --------------------------------------------------
    # 2. Data type specific preprocessing
    # --------------------------------------------------
    if data_type == "rna":
        sc.tl.pca(adata, n_comps=n_pcs, svd_solver="arpack")

    elif data_type == "atac":
        print(f"Processing ATAC data with {adata.n_vars} peaks")

        if atac_n_top_features is not None:
            sc.pp.highly_variable_genes(
                adata,
                n_top_genes=atac_n_top_features,
                flavor="seurat_v3"
            )
            adata = adata[:, adata.var.highly_variable]
            print(f"Selected {atac_n_top_features} top features")

        from sklearn.feature_extraction.text import TfidfTransformer
        from scipy import sparse

        if not sparse.issparse(adata.X):
            adata.X = sparse.csr_matrix(adata.X)

        tfidf = TfidfTransformer(norm='l2')
        adata.X = tfidf.fit_transform(adata.X)

        sc.pp.pca(
            adata,
            n_comps=min(n_pcs, adata.n_vars - 1),
            svd_solver="arpack",
            zero_center=False,
            random_state=42
        )

        adata.obsm['X_lsi'] = adata.obsm['X_pca'].copy()
        adata.uns['pca'] = adata.uns['pca'].copy()
        adata.uns['pca']['params']['zero_center'] = False

    else:
        raise ValueError(f"Unsupported data_type: {data_type}. Use 'rna' or 'atac'")

    # --------------------------------------------------
    # 3. KNN graph construction
    # --------------------------------------------------
    use_rep = 'X_lsi' if data_type == 'atac' else 'X_pca'

    sc.pp.neighbors(
        adata,
        n_neighbors=n_neighbors,
        n_pcs=min(n_pcs, adata.n_vars - 1),
        use_rep=use_rep,
        metric='euclidean' if data_type == 'rna' else 'cosine',
        random_state=42
    )

    # --------------------------------------------------
    # 4. Leiden clustering
    # --------------------------------------------------
    sc.tl.leiden(
        adata,
        resolution=resolution,
        key_added="pseudo_label",
        random_state=42
    )

    if 'true_label' in adata.obs:
        try:
            from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score

            true_labels = adata.obs['true_label'].astype(str).values
            pseudo_labels = adata.obs['pseudo_label'].astype(str).values

            ari = adjusted_rand_score(true_labels, pseudo_labels)
            ami = adjusted_mutual_info_score(true_labels, pseudo_labels)

            print(f"Clustering quality (vs true labels):")
            print(f"  Adjusted Rand Index (ARI): {ari:.4f}")
            print(f"  Adjusted Mutual Info (AMI): {ami:.4f}")

            adata.uns['clustering_metrics'] = {
                'ARI': ari,
                'AMI': ami,
                'n_clusters': len(adata.obs['pseudo_label'].unique()),
                'resolution': resolution
            }

        except Exception as e:
            print(f"Could not compute clustering metrics: {e}")

    from sklearn.preprocessing import LabelEncoder
    le = LabelEncoder()
    pseudo_labels_numeric = le.fit_transform(adata.obs["pseudo_label"])

    adata.uns['pseudo_label_encoder'] = {
        'classes': le.classes_.tolist(),
        'mapping': dict(zip(le.classes_, le.transform(le.classes_)))
    }

    pseudo_labels_tensor = torch.tensor(
        pseudo_labels_numeric,
        dtype=torch.long,
        device=device
    )

    print(f"Generated {len(np.unique(pseudo_labels_numeric))} pseudo clusters")
    print(f"Pseudo label distribution: {np.bincount(pseudo_labels_numeric)}")

    return pseudo_labels_tensor, adata


def leiden_pseudo_label(
    adata,
    use_rep="X_scVI",
    n_neighbors=15,
    resolution=0.5,
    min_cluster_cells=20,
    random_state=42,
    n_id=None,
    key_added="pseudo_label",
):
    """
    Perform Leiden clustering on precomputed latent representation.
    """

    assert use_rep in adata.obsm, f"{use_rep} not found in adata.obsm"

    # --------------------------------------------------
    # 1. KNN graph
    # --------------------------------------------------
    sc.pp.neighbors(
        adata,
        use_rep=use_rep,
        n_neighbors=n_neighbors,
        random_state=random_state,
    )

    # --------------------------------------------------
    # 2. Leiden clustering
    # --------------------------------------------------
    raw_key = f"{key_added}_raw"

    sc.tl.leiden(
        adata,
        resolution=resolution,
        key_added=raw_key,
        random_state=random_state,
    )

    # --------------------------------------------------
    # 3. Small cluster filtering
    # --------------------------------------------------
    labels = adata.obs[raw_key].astype(int).values
    counts = np.bincount(labels)

    filtered = labels.copy()
    for c, cnt in enumerate(counts):
        if cnt < min_cluster_cells:
            filtered[labels == c] = -1

    adata.obs[key_added] = filtered.astype(str)
    adata.obs[key_added] = adata.obs[key_added].astype("category")

    # --------------------------------------------------
    # 4. Stats
    # --------------------------------------------------
    n_clusters = len(set(filtered)) - (1 if -1 in filtered else 0)
    noise_ratio = np.mean(filtered == -1)

    adata.uns[f"{key_added}_info"] = {
        "method": "Leiden",
        "use_rep": use_rep,
        "n_neighbors": n_neighbors,
        "resolution": resolution,
        "min_cluster_cells": min_cluster_cells,
        "n_clusters": n_clusters,
        "noise_ratio": float(noise_ratio),
    }

    print(
        f"[Leiden] clusters={n_clusters}, noise={noise_ratio:.3f}, "
        f"res={resolution}, k={n_neighbors}"
    )

    adata = shift_ood_pseudo_labels(
        adata,
        pseudo_key="pseudo_label",
        n_id=n_id,
    )

    return adata, filtered


def shift_ood_pseudo_labels(adata, pseudo_key="pseudo_label", n_id=None):
    """
    Shift OOD pseudo labels by n_id and convert to categorical.

    Parameters
    ----------
    adata : AnnData
        OOD AnnData
    pseudo_key : str
        Key of pseudo labels in adata.obs
    n_id : int
        Number of ID classes

    Returns
    -------
    adata : AnnData
        Modified AnnData (in-place)
    """

    if n_id is None:
        raise ValueError("n_id must be provided for OOD pseudo labels")

    labels = adata.obs[pseudo_key]

    if hasattr(labels, "cat"):
        labels = labels.cat.codes

    labels = labels.astype(int)

    labels_shifted = labels + n_id

    adata.obs[pseudo_key] = labels_shifted.astype("category")

    return adata


def filter_pseudo_data(adata, dim=100, draw_fig=True):
    labels = adata.obs["pseudo_label"].astype(int).to_numpy()
    Z = adata.obsm["X_scVI"][:, :100]  # encoder latent

    dist_to_center = np.zeros(len(labels))

    for c in np.unique(labels):
        idx = labels == c
        Zc = Z[idx]

        center = Zc.mean(axis=0)
        dist = np.linalg.norm(Zc - center, axis=1)

        dist_to_center[idx] = dist

    keep_mask = np.ones(len(labels), dtype=bool)

    for c in np.unique(labels):
        idx = labels == c
        dist_c = dist_to_center[idx]

        thresh = np.percentile(dist_c, 50)
        keep_mask[idx] = dist_c <= thresh

    adata_filtered = adata[keep_mask].copy()
    adata.obs["dist_to_center"] = dist_to_center

    if draw_fig:
        sc.pl.umap(
            adata,
            color="dist_to_center",
            cmap="viridis",
            frameon=False,
        )
    return adata, adata_filtered


def draw_true_pseudo(adata_filtered):
    if 'true_label' in adata_filtered.obs:
        sc.pl.umap(
            adata_filtered,
            color=["true_label", "pseudo_label"],
            frameon=False,
            ncols=2,
            wspace=0.5,
            title=["True Labels", "Pseudo Labels"]
        )


def clustering_metrics_from_adata(
    adata,
    true_key="cell_type",
    pred_key="pseudo_label",
    ignore_label="-1",
):
    """
    Compute ARI / AMI / NMI / Homogeneity between true labels and pseudo labels.

    Returns
    -------
    metrics : dict
    """

    assert true_key in adata.obs, f"{true_key} not found in adata.obs"
    assert pred_key in adata.obs, f"{pred_key} not found in adata.obs"

    y_true = adata.obs[true_key].astype(str).values
    y_pred = adata.obs[pred_key].astype(str).values

    # --------------------------------------------------
    # 1. Filter noise labels
    # --------------------------------------------------
    mask = y_pred != ignore_label

    if mask.sum() == 0:
        raise ValueError("All pseudo labels are ignored (-1).")

    y_true_f = y_true[mask]
    y_pred_f = y_pred[mask]

    # --------------------------------------------------
    # 2. Compute metrics
    # --------------------------------------------------
    ari = adjusted_rand_score(y_true_f, y_pred_f)
    ami = adjusted_mutual_info_score(y_true_f, y_pred_f)
    nmi = normalized_mutual_info_score(y_true_f, y_pred_f)
    hom = homogeneity_score(y_true_f, y_pred_f)

    metrics = {
        "ARI": ari,
        "AMI": ami,
        "NMI": nmi,
        "HOM": hom,
        "n_cells_used": int(mask.sum()),
        "noise_ratio": float(1 - mask.mean()),
    }

    return metrics


def dataloader_to_adata(
        dataloader,
        label_key="label",  # cell_type
        dtype=np.float32,
):
    """
    Convert a dataloader yielding (x_raw, x, y) into AnnData.

    Parameters
    ----------
    dataloader : torch.utils.data.DataLoader
        Yields (x_raw, x, y)
    label_key : str
        Column name in adata.obs for labels
    dtype : numpy dtype
        Data type for X

    Returns
    -------
    adata : anndata.AnnData
        adata.X is sparse (CSR), adata.obs[label_key] contains y
    """

    process_x_list = []
    X_list = []
    y_list = []

    for batch in dataloader:
        x_raw, x, y = batch

        # ---- x_raw ----
        if torch.is_tensor(x_raw):
            x_raw = x_raw.detach().cpu().numpy()
        X_list.append(x_raw.astype(dtype))

        # ---- x ----
        if torch.is_tensor(x):
            x = x.detach().cpu().numpy()
        process_x_list.append(x.astype(dtype))

        # ---- y ----
        if torch.is_tensor(y):
            y = y.detach().cpu().numpy()
        y_list.append(y)

    # concat
    X = np.concatenate(X_list, axis=0)  # (N, D), dense
    y = np.concatenate(y_list, axis=0)  # (N,)
    process_x = np.concatenate(process_x_list, axis=0)

    # dense -> sparse
    X_sparse = csr_matrix(X)

    # build AnnData
    adata = ad.AnnData(X=X_sparse)
    adata.obs[label_key] = y

    adata.obs["cell_type"] = (
        adata.obs["cell_type"].astype("category")
    )

    adata.obsm['processed_x'] = process_x
    adata.obsm['raw_x'] = X

    return adata


def rna_adata_to_dataloader_with_pseudo(
    adata,
    pseudo_key="pseudo_label",
    label_key="cell_type",
    batch_size=128,
    shuffle=True,
    num_workers=0,
):
    """
    Convert AnnData to DataLoader returning:
    (pseudo_label, data, true_label)

    Data preprocessing:
        normalize_total -> log1p -> scale
    """

    adata = adata.copy()

    adata.X = adata.X.toarray()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    sc.pp.scale(adata)

    X = adata.X

    X = torch.tensor(X, dtype=torch.float32)
    pseudo = adata.obs[pseudo_key]
    true = adata.obs[label_key]

    pseudo = torch.tensor(pseudo.values, dtype=torch.long)
    true = torch.tensor(true.values, dtype=torch.long)

    dataset = TensorDataset(true, X, pseudo)

    loader = DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=shuffle,
        drop_last=False,
    )

    return loader


def atac_adata_to_dataloader(
    adata,
    pseudo_key="pseudo_label",
    label_key="cell_type",
    batch_size=128,
    shuffle=True,
    num_workers=0,
):
    """
    Convert AnnData to DataLoader returning:
    (pseudo_label, data, true_label)

    Data preprocessing:
        normalize_total -> log1p -> scale
    """

    adata = adata.copy()

    adata.X = adata.X.toarray()

    tfidf = TfidfTransformer()
    adata.X = tfidf.fit_transform(adata.X).toarray()
    sc.pp.scale(adata)

    X = adata.X

    X = torch.tensor(X, dtype=torch.float32)

    pseudo = adata.obs[label_key]
    true = adata.obs[label_key]

    pseudo = torch.tensor(pseudo.values, dtype=torch.long)
    true = torch.tensor(true.values, dtype=torch.long)

    dataset = TensorDataset(true, X, pseudo)

    loader = DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=shuffle,
        drop_last=False,
    )

    return loader


def get_pseudo_label(filtered_omics1_ood_loader):
    adata_rna_ood = dataloader_to_adata(filtered_omics1_ood_loader, label_key="cell_type")
    # print(adata_rna_ood.obs['cell_type'].value_counts())

    sc.pp.normalize_total(adata_rna_ood, target_sum=1e4)
    sc.pp.log1p(adata_rna_ood)
    sc.pp.highly_variable_genes(adata_rna_ood, n_top_genes=2000)
    adata_rna_ood = adata_rna_ood[:, adata_rna_ood.var.highly_variable]
    sc.pp.scale(adata_rna_ood, max_value=10)

    sc.tl.pca(adata_rna_ood, svd_solver='arpack')
    sc.pl.pca_variance_ratio(adata_rna_ood, log=True, show=False)  # True
    plt.close()

    sc.pp.neighbors(adata_rna_ood, n_neighbors=50, n_pcs=30)
    sc.tl.umap(adata_rna_ood)

    sc.tl.leiden(adata_rna_ood, resolution=0.1)
    # sc.pl.umap(adata_rna_ood, color=['leiden', 'cell_type'])

    return adata_rna_ood


def get_augmented_data(adata_rna_ood, class_number, omics1_train_loader):
    # print(adata_rna_ood)
    raw_x = adata_rna_ood.obsm['raw_x']
    x = adata_rna_ood.obsm['processed_x']
    # print(x.shape)
    y = adata_rna_ood.obs['leiden']
    y = y.astype('int') + class_number

    raw_x = torch.tensor(raw_x)
    x = torch.tensor(x)
    y = torch.tensor(y)

    filtered_ood_dataset = BuildDataset(raw_x, x, y)

    train_id_dataset = omics1_train_loader.dataset

    from torch.utils.data import ConcatDataset

    merged_dataset = ConcatDataset([train_id_dataset, filtered_ood_dataset])

    merged_loader = DataLoader(
        dataset=merged_dataset, batch_size=128, shuffle=True, drop_last=True
    )

    return merged_loader


def get_class_number_from_loader(loader):
    labels = []
    for _, _, y in loader:
        y = y.view(-1)
        labels.append(y)
    labels = torch.cat(labels)

    return int(labels.max().item() + 1)


def expansion_fine_tuning(E_rna, D_rna, Classifier, optimizer_E_rna, optimizer_D_rna, scheduler_E_rna, scheduler_D_rna,
                           merged_loader, omics1_test_id_loader, num_epochs, device):
    num_classes_new = get_class_number_from_loader(merged_loader)

    Classifier_expanded = expand_classifier(Classifier, num_classes_new)

    optimizer_fc = optim.Adam(Classifier_expanded.parameters(), lr=0.001, weight_decay=1e-4)
    scheduler_fc = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer_fc, T_max=100)

    E_rna, D_rna, Classifier_expanded = train_single_omics(E_rna, D_rna, Classifier_expanded,
                                                                                  optimizer_E_rna, optimizer_D_rna,
                                                                                  optimizer_fc, scheduler_E_rna,
                                                                                  scheduler_D_rna, scheduler_fc,
                                                                                  merged_loader, omics1_test_id_loader,
                                                                                  num_epochs=num_epochs, device=device)

    return E_rna, D_rna, Classifier_expanded


def get_predict_for_evaluation_with_single_omics(E_rna, Classifier_expanded, omics1_mixed_loader, device):
    E_rna.eval()
    Classifier_expanded.eval()

    correct = 0
    total = 0
    y_true = []
    y_pred = []
    X_raw = []
    emb = []
    scores = []

    with torch.no_grad():
        for x_raw, x_rna, y_rna in omics1_mixed_loader:
            x_rna = x_rna.to(device)
            y_rna = y_rna.to(device)

            z_rna, _, _ = E_rna(x_rna)

            logits = Classifier_expanded(z_rna)

            smax = to_np(F.softmax(logits, dim=1))

            score = np.max(smax, axis=1)  # np.max(smax, axis=1)

            _, predicted = torch.max(logits.data, 1)

            emb.append(z_rna.detach().cpu())
            y_true.append(y_rna.detach().cpu())
            y_pred.append(predicted.detach().cpu())
            X_raw.append(x_rna.detach().cpu())
            scores.append(score[:len(x_rna)])

        emb = torch.cat(emb, dim=0)  # [:max_samples]
        y_true = torch.cat(y_true, dim=0)
        y_pred = torch.cat(y_pred, dim=0)
        X_raw = torch.cat(X_raw, dim=0)
        Score = np.concatenate(scores, axis=0)
    # print(f'Accuracy of the model on the test images: {accuracy:.2f}%')
    return emb.numpy(), y_true.numpy(), y_pred.numpy(), X_raw.numpy(), Score


def get_expanded_prediction_results(E_rna, Classifier_expanded, omics1_mixed_loader, device, gene_name, class_number):
    emb, y_true, y_pred, X_raw, ood_score = get_predict_for_evaluation_with_single_omics(E_rna, Classifier_expanded,
                                                                                         omics1_mixed_loader, device)

    adata = sc.AnnData(X=X_raw)

    adata.var["gene_symbol"] = gene_name

    adata.obs['y_true'] = pd.Series(y_true, index=adata.obs_names).astype(str)
    adata.obs['y_pred'] = pd.Series(y_pred, index=adata.obs_names).astype(str)
    adata.obs['predicted_score'] = pd.Series(ood_score, index=adata.obs_names).astype(float)

    adata.obsm['X_scDiscovery'] = emb

    y_true_new = adata.obs['y_true'].copy()
    y_pred_new = adata.obs['y_pred'].copy()

    known_classes = [str(i) for i in range(class_number)]  # ['0','1','2','3','4','5','6']

    y_true_new = y_true_new.where(y_true_new.isin(known_classes), 'Unknown')
    y_pred_new = y_pred_new.where(y_pred_new.isin(known_classes), 'Unknown')

    adata.obs['y_true_with_unknown'] = y_true_new.astype(str)
    adata.obs['y_pred_with_unknown'] = y_pred_new.astype(str)

    return adata


def dynamic_novel_cell_type_expansion(E_rna, D_rna, Classifier, optimizer_E_rna, optimizer_D_rna, scheduler_E_rna, scheduler_D_rna,
                                                              filtered_omics1_ood_loader, omics1_train_loader, omics1_test_id_loader, omics1_mixed_loader, gene_name, class_number, dataset_name, num_epochs, device):


    adata_rna_ood = get_pseudo_label(filtered_omics1_ood_loader)
    merged_loader = get_augmented_data(adata_rna_ood, class_number, omics1_train_loader)

    E_rna, D_rna, Classifier_expanded = expansion_fine_tuning(E_rna, D_rna, Classifier, optimizer_E_rna,
                                                              optimizer_D_rna, scheduler_E_rna, scheduler_D_rna,
                                                              merged_loader, omics1_test_id_loader, num_epochs, device)

    adata = get_expanded_prediction_results(E_rna, Classifier_expanded, omics1_mixed_loader, device, gene_name,
                                            class_number)

    adata.write(f'./output/scDiscovery_{dataset_name}.h5ad')

    return adata, E_rna, D_rna, Classifier_expanded


def cell_type_mapping(
    adata,
    global_categories_id,
    mapping_ood,
    y_true_col="y_true",
    y_pred_col="y_pred",
    true_name_col="y_true_name",
    pred_name_col="y_pred_name",
    novel_prefix="Pred Novel Type",
    inplace=True,
    verbose=False
):
    import pandas as pd

    if not inplace:
        adata = adata.copy()

    id_to_name = {
        int(i): name
        for i, name in enumerate(global_categories_id)
    }

    known_ids = set(id_to_name.keys())

    true_mapping = {
        **id_to_name,
        **{int(k): v for k, v in mapping_ood.items()}
    }

    y_pred_int = pd.to_numeric(
        adata.obs[y_pred_col],
        errors="coerce"
    ).astype("Int64")

    y_true_int = pd.to_numeric(
        adata.obs[y_true_col],
        errors="coerce"
    ).astype("Int64")

    pred_ids = y_pred_int.dropna().astype(int).unique()

    ood_pred_ids = sorted([
        int(x) for x in pred_ids
        if int(x) not in known_ids
    ])

    pred_ood_mapping = {
        pred_id: f"{novel_prefix} {i}"
        for i, pred_id in enumerate(ood_pred_ids, start=1)
    }

    adata.obs[true_name_col] = y_true_int.map(true_mapping)

    adata.obs[pred_name_col] = y_pred_int.map(id_to_name)

    ood_mask = y_pred_int.isin(pred_ood_mapping.keys())
    adata.obs.loc[ood_mask, pred_name_col] = y_pred_int[ood_mask].map(pred_ood_mapping)

    adata.obs[true_name_col] = adata.obs[true_name_col].fillna(
        y_true_int.map(lambda x: f"unknown true id {x}" if pd.notna(x) else pd.NA)
    )

    adata.obs[pred_name_col] = adata.obs[pred_name_col].fillna(
        y_pred_int.map(lambda x: f"unknown pred id {x}" if pd.notna(x) else pd.NA)
    )

    if verbose:
        print("id_to_name:")
        print(id_to_name)

        print("\ntrue_mapping:")
        print(true_mapping)

        print("\npred_ood_mapping:")
        print(pred_ood_mapping)

        print("\ncheck:")
        print(
            adata.obs[
                [y_true_col, true_name_col, y_pred_col, pred_name_col]
            ].head(10)
        )

    return adata