import torch
import scanpy as sc
import numpy as np
import pandas as pd
from torch.utils.data import TensorDataset, DataLoader, Dataset, Subset, ConcatDataset
from sklearn.feature_extraction.text import TfidfTransformer
from data_configs import *
import warnings
warnings.filterwarnings("ignore")
from logger import *
import os


class BuildDataset(Dataset):
    def __init__(self, data_raw, data, labels):
        self.data_raw = data_raw
        self.data = data
        self.labels = labels

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        return self.data_raw[idx] , self.data[idx], self.labels[idx]


def load_data(dataset_name, setting='setting1', omics=1):
    dataset = dataset_dict[dataset_name]
    modalities = ood_setting_dict[dataset_name]["modalities"]
    cell_type_key = ood_setting_dict[dataset_name][setting]["cell_type_key"]
    batch_key = ood_setting_dict[dataset_name][setting]["batch_key"]
    remove_type = ood_setting_dict[dataset_name][setting]["remove_type"]

    if omics == 1:
        omics1 = sc.read_h5ad(f"../../data/sc/{dataset}.h5ad")
        omics1.var_names_make_unique()

        if 'RNA' in modalities:
            sc.pp.highly_variable_genes(
                omics1,
                n_top_genes=2000,
                flavor="seurat_v3"
            )
            omics1 = omics1[:, omics1.var.highly_variable].copy()

        omics2 = None

        if dataset_name == 'peng':
            gene_name = omics1.var_names.astype(str)
        elif dataset_name == 'mice_zebrafish':
            omics1.var['ensembl_id'] = omics1.var_names.astype(str)
            omics1.var_names = omics1.var['feature_name'].astype(str)
            omics1.var_names_make_unique()
            gene_name = omics1.var_names.astype(str)
        elif dataset_name == 'Zeisel_p20_p29':
            gene_name = omics1.var_names.astype(str)
        else:
            gene_name = omics1.var_names.astype(str)
        omics1_train_loader, omics1_test_id_loader, omics1_test_ood_loader, global_categories_id, global_categories_ood, mapping_ood = split_train_test_by_batch(
            omics1, omics2, dataset_name, setting, cell_type_key, batch_key, remove_type, modalities)
        return omics1_train_loader, omics1_test_id_loader, omics1_test_ood_loader, global_categories_id, global_categories_ood, mapping_ood, gene_name
    elif omics == 2:
        omics1 = sc.read_h5ad(f"../Data/sc/Data/{dataset}-{modalities[0]}.h5ad")
        omics2 = sc.read_h5ad(f"../Data/sc/Data/ATAC_gene_activity/{dataset}-{modalities[1]}_gene_activity.h5ad")
        omics1.var_names_make_unique()
        omics2.var_names_make_unique()

        omics1_train_loader, omics1_test_id_loader, omics1_test_ood_loader, omics2_train_loader, omics2_test_id_loader, omics2_test_ood_loader, global_categories_id, global_categories_ood, mapping_ood = split_train_test_by_batch(omics1, omics2, dataset_name, setting)

        return omics1_train_loader, omics1_test_id_loader, omics1_test_ood_loader, omics2_train_loader, omics2_test_id_loader, omics2_test_ood_loader, global_categories_id, global_categories_ood, mapping_ood


def split_train_test_by_batch(omics1, omics2, dataset_name, setting, cell_type_key, batch_key, remove_type, modalities):
    # split train data and test data
    # omics1
    ood_type = ood_setting_dict[dataset_name][setting]['ood_type']
    train_batch = ood_setting_dict[dataset_name][setting]['train_batch']
    test_batch = ood_setting_dict[dataset_name][setting]['test_batch']

    split_by_batch_logger = create_logger(name='Split Train and Test Data by Batch', ch=True, fh=False,
                                          levelname=logging.INFO, overwrite=False)
    split_by_batch_logger.info(f'Train: {train_batch}')
    split_by_batch_logger.info(f'Test: {test_batch}')

    mask_id_type = ~omics1.obs[cell_type_key].isin(ood_type)
    mask_train_batch = omics1.obs[batch_key].isin(train_batch)

    omics1_train = omics1[mask_id_type & mask_train_batch].copy()
    omics1_test_id = omics1[~mask_train_batch & mask_id_type].copy()
    omics1_test_ood = omics1[~mask_train_batch & ~mask_id_type].copy()

    # remove 指定类别
    omics1_test_id = omics1_test_id[~omics1_test_id.obs[cell_type_key].isin(remove_type)].copy()
    omics1_train = omics1_train[~omics1_train.obs[cell_type_key].isin(remove_type)].copy()

    # omics2
    if omics2 != None:
        ood_type = ood_setting_dict[dataset_name][setting]['ood_type']
        train_batch = ood_setting_dict[dataset_name][setting]['train_batch']
        test_batch = ood_setting_dict[dataset_name][setting]['test_batch']
        modality = ood_setting_dict[dataset_name]['modalities']

        mask_id_type = ~omics2.obs[cell_type_key].isin(ood_type)
        mask_train_batch = omics2.obs[batch_key].isin(train_batch)

        omics2_train = omics2[mask_id_type & mask_train_batch].copy()
        # omics2_test = omics2[~mask_train_batch].copy()
        omics2_test_id = omics2[~mask_train_batch & mask_id_type].copy()
        omics2_test_ood = omics2[~mask_train_batch & ~mask_id_type].copy()

        # remove 指定类别
        omics2_test_id = omics2_test_id[~omics2_test_id.obs[cell_type_key].isin(remove_type)].copy()

        global_categories_id = build_global_categories(
            omics1_train.obs[cell_type_key],
            omics2_train.obs[cell_type_key],
        )

        global_categories_ood = build_global_categories(
            omics1_test_ood.obs[cell_type_key],
            omics2_test_ood.obs[cell_type_key],
        )
    else:
        global_categories_id = build_global_categories(
            omics1_train.obs[cell_type_key],
        )

        global_categories_ood = build_global_categories(
            omics1_test_ood.obs[cell_type_key],
        )

        print(global_categories_id)

    # Make Datasets
    # RNA
    omics1_train_loader, mapping_id = make_dataloader(omics1_train, global_categories_id, process_adata, cell_type_key, modalities, is_ood=False, n_id=len(global_categories_id))
    if omics1_test_id.shape[0] != 0:
        omics1_test_id_loader, _ = make_dataloader(omics1_test_id, global_categories_id, process_adata, cell_type_key, modalities, is_ood=False, n_id=len(global_categories_id))
    else:
        omics1_test_id_loader = None
    omics1_test_ood_loader, mapping_ood = make_dataloader(omics1_test_ood, global_categories_ood, process_adata, cell_type_key, modalities, is_ood=True, n_id=len(global_categories_id))

    if omics2 == None:
        return omics1_train_loader, omics1_test_id_loader, omics1_test_ood_loader, global_categories_id, global_categories_ood, mapping_ood
    else:
        # ATAC
        omics2_train_loader, _ = make_dataloader(omics2_train, global_categories_id, process_atac_adata, cell_type_key, modalities, is_ood=False, n_id=len(global_categories_id))
        omics2_test_id_loader, _ = make_dataloader(omics2_test_id, global_categories_id, process_atac_adata, cell_type_key, modalities, is_ood=False, n_id=len(global_categories_id))
        omics2_test_ood_loader, _ = make_dataloader(omics2_test_ood, global_categories_ood, process_atac_adata, cell_type_key, modalities, is_ood=True, n_id=len(global_categories_id))

        return omics1_train_loader, omics1_test_id_loader, omics1_test_ood_loader, omics2_train_loader, omics2_test_id_loader, omics2_test_ood_loader, global_categories_id, global_categories_ood, mapping_ood


def make_dataloader(adata, global_categories, data_process, cell_type_key, modalities, is_ood=False, n_id=None):
    data_raw, data, label, mapping = data_process(
        adata,
        categories=global_categories,
        cell_type_key=cell_type_key,
        modality_key=modalities,
        is_ood=is_ood,
        n_id=n_id
    )
    omics_dataset = BuildDataset(data_raw, data, label)

    omics_loader = DataLoader(
        dataset=omics_dataset, batch_size=128, shuffle=True, drop_last=False
    )

    return omics_loader, mapping


def build_global_categories(*cell_type_series_list):
    all_types = pd.concat(cell_type_series_list)
    categories = pd.Index(sorted(all_types.unique()))
    return categories


# def encode_labels_with_categories(cell_type_series, categories):
#     cat = pd.Categorical(cell_type_series, categories=categories)
#     if (cat.codes < 0).any():
#         raise ValueError("Found unknown cell type not in global categories")
#     return cat.codes.astype(int)


def encode_labels_with_categories(
    cell_type_series,
    categories,
    *,
    is_ood: bool = False,
    n_id: int = None,
    return_mapping: bool = False,
):
    import pandas as pd

    cell_type_series = pd.Series(cell_type_series)

    if not is_ood:
        cat = pd.Categorical(cell_type_series, categories=categories)
        if (cat.codes < 0).any():
            raise ValueError("Found unknown cell type in ID data")

        labels = cat.codes.astype(int)

        if return_mapping:
            mapping = {i: ct for i, ct in enumerate(categories)}
            return labels, mapping

        return labels

    # -------- OOD --------
    if n_id is None:
        raise ValueError("n_id must be provided when is_ood=True")

    ood_types = pd.unique(cell_type_series)

    # forward mapping
    ood_mapping = {
        ct: n_id + i
        for i, ct in enumerate(ood_types)
    }

    labels = cell_type_series.map(ood_mapping).astype(int).values

    if return_mapping:
        # inverse mapping: label -> cell type
        inv_mapping = {v: k for k, v in ood_mapping.items()}
        return labels, inv_mapping

    return labels


def process_adata(adata_rna, categories, cell_type_key, modality_key, is_ood=False, n_id=None):
    #############
    adata_rna.X = adata_rna.X.toarray()
    data_raw = torch.tensor(adata_rna.X).float()
    #############

    if 'ADT' in modality_key:
        sc.pp.normalize_total(adata_rna, target_sum=1e4)
        sc.pp.log1p(adata_rna)

    codes, mapping = encode_labels_with_categories(adata_rna.obs[cell_type_key], categories, is_ood=is_ood, n_id=n_id, return_mapping=True)
    data = torch.tensor(adata_rna.X).float()
    label = torch.tensor(codes).long()
    return data_raw, data, label, mapping


def concat_datasets(dataloader_rna, dataloader_atac):
    dataset_rna = dataloader_rna.dataset
    dataset_atac = dataloader_atac.dataset

    mixed_dataset = ConcatDataset([dataset_rna, dataset_atac])

    mixed_loader = DataLoader(
        dataset=mixed_dataset, batch_size=128, shuffle=True, drop_last=False
    )

    return mixed_loader
