import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score, normalized_mutual_info_score, homogeneity_score
import torch
from scDiscovery_architecture import reparameterize
from sklearn.metrics import (
    classification_report,
    f1_score,
    accuracy_score,
    silhouette_score,
    calinski_harabasz_score
)


def evaluate_discovery_potential(adata, y_true, y_pred, embed_key='X_glue', unknown_label='Unknown'):
    y_true_arr = adata.obs[y_true].astype(str).values
    y_pred_arr = adata.obs[y_pred].astype(str).values
    emb = adata.obsm[embed_key]

    metrics = {}

    metrics['Multi-F1'] = f1_score(y_true_arr, y_pred_arr, average='weighted')
    metrics['Multi-Accuracy'] = accuracy_score(y_true_arr, y_pred_arr)

    true_is_unk = (y_true_arr == unknown_label)
    pred_is_unk = (y_pred_arr == unknown_label)

    metrics['F1'] = f1_score(true_is_unk, pred_is_unk, pos_label=True)

    if pred_is_unk.sum() > 0:
        metrics['Precision'] = (true_is_unk & pred_is_unk).sum() / pred_is_unk.sum()
    else:
        metrics['Precision'] = 0.0

    if true_is_unk.sum() > 0:
        metrics['Recall'] = (true_is_unk & pred_is_unk).sum() / true_is_unk.sum()
    else:
        metrics['Recall'] = 0.0

    return metrics


def _calculate_cluster_metrics(emb, mask, prefix):
    import scanpy as sc
    res = {}

    if mask.sum() < 20:
        res[f'{prefix}_ASW'] = 0.0
        res[f'{prefix}_N_Clusters'] = 0
        return res

    emb_subset = emb[mask]

    try:
        temp_adata = sc.AnnData(X=np.zeros((emb_subset.shape[0], 1)))
        temp_adata.obsm['X_emb'] = emb_subset

        sc.pp.neighbors(temp_adata, use_rep='X_emb', n_neighbors=15)
        sc.tl.leiden(temp_adata, resolution=0.5, key_added='sub')

        labels = temp_adata.obs['sub'].values
        n_clusters = len(np.unique(labels))

        res[f'{prefix}_N_Clusters'] = n_clusters

        if n_clusters > 1:
            res[f'{prefix}_ASW'] = silhouette_score(emb_subset, labels)
            res[f'{prefix}_CH'] = calinski_harabasz_score(emb_subset, labels)
        else:
            res[f'{prefix}_ASW'] = 0.0
            res[f'{prefix}_CH'] = 0.0

    except Exception as e:
        print(f"Clustering failed for {prefix}: {e}")
        res[f'{prefix}_ASW'] = 0.0

    return res


def calculate_discovery_asw(X_emb, n_neighbors=15, resolution=0.5):
    if X_emb.shape[0] < 10:
        return 0.0

    try:
        temp_adata = sc.AnnData(X=np.zeros((X_emb.shape[0], 1)))
        temp_adata.obsm['X_emb'] = X_emb

        sc.pp.neighbors(temp_adata, use_rep='X_emb', n_neighbors=n_neighbors)
        sc.tl.leiden(temp_adata, resolution=resolution, key_added='leiden_temp')

        labels = temp_adata.obs['leiden_temp'].values
        n_clusters = len(np.unique(labels))

        if n_clusters > 1:
            score = silhouette_score(X_emb, labels)
            return score
        else:
            return 0.0

    except Exception as e:
        print(f"{e}")
        return 0


def evaluate_model_on_novel_cell_type(E_rna, D_rna, Classifier, omics1_train_loader, device):
    E_rna.eval()
    Classifier.eval()

    y_true = []
    y_pred = []
    recon_data = []
    emb = []

    correct = 0
    total = 0

    with torch.no_grad():
        for _, x_rna, y_rna in omics1_train_loader:
            x_rna = x_rna.to(device)
            y_rna = y_rna

            z_rna, mu_rna, logvar_rna = E_rna(x_rna)

            ###############################################
            re_z_rna = reparameterize(mu_rna, logvar_rna)

            # ---------------- decoder -----------------
            x_hat_rna = D_rna(re_z_rna)

            ###############################################

            logits = Classifier(z_rna)

            _, predicted = torch.max(logits.data, 1)

            total += logits.size(0)
            correct += (predicted.cpu() == y_rna).sum().item()

            y_true.append(y_rna)
            y_pred.append(predicted.cpu())
            recon_data.append(x_hat_rna.cpu())
            emb.append(z_rna.detach().cpu())

    accuracy = 100 * correct / total

    y_true = torch.cat(y_true, dim=0).numpy()
    y_pred = torch.cat(y_pred, dim=0).numpy()
    recon_data_rna = torch.cat(recon_data, dim=0)
    Emb = torch.cat(emb, dim=0).numpy()

    ari = adjusted_rand_score(y_true, y_pred)
    ami = adjusted_mutual_info_score(y_true, y_pred)
    nmi = normalized_mutual_info_score(y_true, y_pred)
    hom = homogeneity_score(y_true, y_pred)

    metrics = {
        "ARI": ari,
        "AMI": ami,
        "NMI": nmi,
        "HOM": hom
    }

    return metrics, accuracy, recon_data_rna, torch.tensor(y_true), torch.tensor(y_pred), Emb
