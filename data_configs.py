dataset_dict = {
    "PBMC": "10x-Multiome-Pbmc10k",
    "Mouse-Cortex": "Chen-2019",
    "Mouse-Skin": "Ma-2020",
    "Human-Kidney-RNA": "UP_HK-RNA",
    "Human-Kidney-ATAC": "UP_HK-ATAC",
    "Human-Kidney": "UP_HK",
    "Mouse-MOp-RNA": "Yao-2021-RNA",
    "Mouse-MOp-ATAC": "Yao-2021-ATAC",
    "BMMC": "CITE-BMMC",
    "Pancreas": "Pancreas",
    "TMS": "TMS_Pancreas",
    "PancreasHumanMouse": "Pancreas_Human_Mouse_Baron",
    "ImmuneHumanMouse": "Immune_Bone_Marrow_Human_Mouse",
    "Baron_Human_Mouse_Normalized": "Baron_Human_Mouse_Normalized",
    "Lung": "Lung_atlas_public",
    "PurifiedPBMCDataset": "PurifiedPBMCDataset",
    # "TS_Bone_Marrow": "TS_Bone_Marrow",
    # "TS_Skin": "TS_Skin",
    "Dorsal_Midbrain": "Dorsal_midbrain_cell_bin",
    "BGI": "BGI_all_coarse",
    "CL_RNA": "CL_RNA",
    "CL_ATAC": "CL_ATAC",
    "PBMC_Time": "pbmc_seurat_v4",
    "PBMC_Sub": "pbmc_seurat_v4",
    "Immune_PBMC_Human_Mouse": "Immune_PBMC_Human_Mouse",
    # "Colorectal_Cancer": "external_lee_natgenet_2020_32451460",
    # "Breast_Cancer": "external_wu_natgenet_2021_34493872",
    # "RCC": "external_zhang_procnatlacadsciusa_2021_34099557",
    "HGSOC": "external_slyper_natmed_2020_32405060",
    "BMMC_Kidney": "BMMC_Kidney",
    "TMS_Fat_Limb": "TMS-Fat-Limb",
    "TMS_marrow_spleen": "TMS-marrow-spleen",
    "Celegans": "celegans",
    "Mus_Kidney": "ALIGNED_Mus_musculus_Kidney",
    "Giraddi_10x": "Giraddi_10x",
    "Han": "Han",
    "Hochane": "Hochane",
    "mouse_heart_pancreas": "mouse_heart_pancreas_with_disease",
    "mice_zebrafish": "mouse_heart_pancreas",  # mice_zebrafish
    "mice_zebrafish_nonCRRL": "mouse_heart_pancreas",  # mice_zebrafish
    "mice_zebrafish_nonOE": "mouse_heart_pancreas",  # mice_zebrafish
    "human_cortex_EaFet_Adol": "human_cortex_EaFet_Adol",
    "Zeisel_p20_p23": "Zeisel_2018_p20_p23",
    "Zeisel_p20_p29": "Zeisel_2018_p20_p29",
    "Zeisel_p20_p29_nonCRRL": "Zeisel_2018_p20_p29",
    "Zeisel_p20_p29_nonOE": "Zeisel_2018_p20_p29",
    "Retina_Human_Mouse": "retina_human_mouse",
    "Retina_Human_Mouse_v2": "Retina_Human_Mouse_Lyu_Macosko",
    "Retina_Human_Mouse_v3": "Retina_Human_Mouse_Lukowski_Shekhar",
    "Thymus_Cao_ALIGNED": "Thymus_Human_Mouse_Cao_ALIGNED",
    "Liver_Human_Monkey": "liver_human_monkey_merged",
    "Liver_Human_Zebrafish": "liver_human_zebrafish_merged",
    "Testis_Human_Platypus": "testis_human_platypus_merged",
    "Brain_Human_Turtle": "brain_human_turtle_merged",
    "Brain_Human_Mouse": "brain_human_mouse_merged",
    "Liver_Human_Mouse": "liver_human_mouse_merged",
    "Glia_Human_Rhesus": "Glia_human_rhesus_merged",
    "Mouse_Pancreas_Subtissue": "tabula-muris-senis-droplet-processed-official-annotations-Pancreas",
    "ALIGNED_Mus_musculus_Mammary_Gland": "ALIGNED_Mus_musculus_Mammary_Gland",
    "TS_Prostate": "TS_Trachea",
    "TMS_Thymus": "tabula-muris-senis-droplet-processed-official-annotations-Thymus",
    "TS_Pancreas": "TS_Pancreas",
    # "TMS-facs-SCAT": "TMS-facs-SCAT",
    "TMS-facs-Trachea": "TMS-facs-Trachea",
    # "TMS-facs-MAT": "TMS-facs-MAT",
    "TMS-droplet-Pancreas": "TMS-droplet-Pancreas",
    # "TMS-facs-Bladder": "",
    "TMS-facs-Limb_Muscle": "TMS-facs-Limb_Muscle",
    "TMS-facs-Lung": "TMS-facs-Lung",
    "TMS-facs-Mammary_Gland": "TMS-facs-Mammary_Gland",
    "TMS-facs-Spleen": "TMS-facs-Spleen",
    # "qian": "external_qian_cellres_2020_32561858",
    "peng": "external_peng_cellres_2019_31273297",
    # "nath": "external_nath_natcommun_2021_34031395",
    "bi": "external_bi_cancercell_2021_33711272",
    "bassez": "external_bassez_natmed_2021_33958794",
    # "Kidney_HeartA": "Kidney_HeartA",
    # "Kidney_Liver": "Kidney_Liver",
    "BoneMarrowA_BoneMarrowB": "BoneMarrowA_BoneMarrowB",
    "MosA1_MosA2": "MosA1_MosA2",
    # "LargeIntestineA_LargeIntestineB": "LargeIntestineA_LargeIntestineB",
    # "LungA_LungB": "LungA_LungB",
    # "MosM1_MosM2": "MosM1_MosM2",
    # "MosP1_MosP2": "MosP1_MosP2",
    # "WholeBrainA_WholeBrainB": "WholeBrainA_WholeBrainB",
    "TMS-facs-GAT": "TMS-facs-GAT",
    # "TMS-facs-Liver": "TMS-facs-Liver",
    # "TMS-droplet-Large_Intestine": "TMS-droplet-Large_Intestine",
    # "TMS-facs-BAT": "TMS-facs-BAT",
    # "TMS-facs-Diaphragm": "TMS-facs-Diaphragm",
    "TMS-facs-Kidney": "TMS-facs-Kidney",
    "TMS-facs-Aorta": "TMS-facs-Aorta",
    "Cao_2020_Eye": "Cao_2020_Eye",
    "Hu": "Hu",
    # "TMS-droplet-Kidney": "TMS-droplet-Kidney",
    # "TMS-Skin": "TMS-droplet-Skin",
    # "TS_Bone_Marrow": "TS_Bone_Marrow",
    "TS_Skin": "TS_Skin",
    "cd34_multiome_rna": "cd34_multiome_rna",
    "bm_multiome_rna": "bm_multiome_rna",
    "TS_Thymus": "TS_Thymus",
    "PBMC_ADT": "pbmc_seurat_v4_adt",
    "BMMC_RNA": "CITE_BMMC_RNA",
    "BMMC_ADT": "CITE_BMMC_ADT",
    "Ma_RNA": "Ma-2020-RNA",
    "Ma_ATAC": "Ma-2020-ATAC_gene_activity",
    "Chen_RNA": "Chen-2019-RNA",
    "Chen_ATAC": "Chen-2019-ATAC",
    "PBMC_10x_RNA": "10x-Multiome-Pbmc10k-RNA",
    "PBMC_10x_ATAC": "10x-Multiome-Pbmc10k-ATAC",
    "Cao_2020_Pancreas": "Cao_2020_Pancreas",
    "Alemany_Fin": "Alemany_Fin",
    "mice_kidney_heart": "mouse_heart_kidney",
    "mice_liver_heart": "mouse_heart_liver",
    "mice_blood_heart": "mouse_heart_blood",
    "PBMC_ATAC_ASAP_seq": "PBMC_ATAC_ASAP_seq",
    "MouseKidney_ATAC": "MouseKidney_ATAC",
    "Human_Bone_Marrow_CITE_ADT": "Human_Bone_Marrow_CITE_ADT",
    "hnscc_rna": "Roider_et_al_BNHL_panel1_FL_rLN",
    "HNSCC_RNA": "hnscc_rna_test_split_rawX",  # hnscc_rna_test_split_rawX hnscc_rna_filtered_2_rawX
    "HNSCC_ADT": "hnscc_adt_test_split",  # HNSCC-ADT  hnscc_adt_test_split  hnscc_adt_filtered_2
    "HNSCC": "HNSCC"
}
ood_setting_dict = {
    "PBMC": {
        "modalities": ["RNA", "ATAC"],
        "batches": ['e0def004-9e30-4a3b-9a65-007110f3a1f2', 'f6c0f811-2fb8-4989-b796-37c14b055517', '5028f75a-8c09-4155-a232-ad7dbfa6042e', '8c570254-4bef-48d8-bd79-c812f60835a5', '8213a3f7-2437-4e8a-b836-caec33df901d'],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_type",
            "batch_key": "batch",
            "remove_type": [],
            "ood_type": ["Memory B", "Naive B", "cDC", "gdT"],
            "split_type": "batch", # "cell_type" or "batch"
            "train_batch": ['e0def004-9e30-4a3b-9a65-007110f3a1f2', 'f6c0f811-2fb8-4989-b796-37c14b055517', '5028f75a-8c09-4155-a232-ad7dbfa6042e', '8c570254-4bef-48d8-bd79-c812f60835a5', '8213a3f7-2437-4e8a-b836-caec33df901d'],
            "test_batch": ['8c570254-4bef-48d8-bd79-c812f60835a5', '8213a3f7-2437-4e8a-b836-caec33df901d']}
    },
    "Mouse-Skin": {
        "modalities": ["RNA", "ATAC_gene_activity"],
        "batches": ['e0def004-9e30-4a3b-9a65-007110f3a1f2', 'f6c0f811-2fb8-4989-b796-37c14b055517', '5028f75a-8c09-4155-a232-ad7dbfa6042e', '8c570254-4bef-48d8-bd79-c812f60835a5', '8213a3f7-2437-4e8a-b836-caec33df901d'],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_type",
            "batch_key": "batch",
            "ood_type": ['Infundibulum', 'Spinous', 'TAC-1', 'TAC-2'],  #
            "remove_type": [],  #
            "train_batch": ['53', '54'],  #
            "test_batch": ['55', '56'],  #
        }
    },
    "HNSCC": {
        "modalities": ["RNA", "ADT"],
        "batches": [],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_type",
            "batch_key": "patient",  # patient_tissue
            "ood_type": ['Macrophage', 'CD8_Exhausted'],  # Macrophage Mast cell
            "remove_type": [],
            "train_batch": ['indiv_2'],  # 128_tumor
            "test_batch": ['indiv_1'],  # 129_tumor
        },
    },
    "HNSCC_RNA": {
        "modalities": ["RNA"],
        "batches": [],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_type",
            "batch_key": "split",
            "ood_type": ['Macrophage', 'CD8_Exhausted'],
            "remove_type": [],
            "train_batch": ['train'],
            "test_batch": ['test'],
        },
        "setting2": {
            "info": "",
            "cell_type_key": "cell_type",
            "batch_key": "patient",  # patient_tissue
            "ood_type": ['Macrophage', 'CD8_Exhausted'],  # Macrophage Mast cell
            "remove_type": [],
            "train_batch": ['indiv_2'],  # 128_tumor
            "test_batch": ['indiv_1'],  # 129_tumor
        },
        "setting3": {
            "info": "",
            "cell_type_key": "cell_type",
            "batch_key": "patient_tissue",
            "ood_type": ['Mast cell', 'CD8_Exhausted'],  # 'B cell'  'CD8 T cell',
            "remove_type": [],  # , 'Endothelial cell', 'Fibroblast', 'Salivary gland' 'pDC'
            "train_batch": ['129_lymph node'],  # , '5' , '128_tumor'  '128_lymph node',
            "test_batch": ['129_tumor'],  #
        },
    },
    "HNSCC_ADT": {
        "modalities": ["ADT"],
        "batches": [],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_type",
            "batch_key": "split",
            "ood_type": ['Macrophage', 'CD8_Exhausted'],  # 'B cell'  'CD8 T cell',
            "remove_type": [],  # , 'Endothelial cell', 'Fibroblast', 'Salivary gland' 'pDC'
            "train_batch": ['train'],  # , '5'
            "test_batch": ['test'],  #
        },
        "setting2": {
            "info": "",
            "cell_type_key": "cell_type",
            "batch_key": "patient",  # patient_tissue
            "ood_type": ['Macrophage', 'CD8_Exhausted'],  # Macrophage  Mast cell
            "remove_type": [],
            "train_batch": ['indiv_2'],  # 128_tumor
            "test_batch": ['indiv_1'],  # 129_tumor
        },
        "setting3": {
            "info": "",
            "cell_type_key": "cell_type",
            "batch_key": "patient_tissue",
            "ood_type": ['Mast cell', 'CD8_Exhausted'],  # 'B cell'  'CD8 T cell',
            "remove_type": [],  # , 'Endothelial cell', 'Fibroblast', 'Salivary gland' 'pDC'
            "train_batch": ['129_lymph node'],  # , '5' , '128_tumor'  '128_lymph node',
            "test_batch": ['129_tumor'],  #
        },
    },
    "Mouse-Cortex": {
        "modalities": ["RNA", "ATAC"],
        "setting1": ['E5Tshz2', 'Endo', 'InN', 'OPC', 'Mic', 'OliI', 'Peri'],
    },
    # "Mouse-Skin": {
    #     "modalities": ["RNA", "ATAC"],
    #     "batches": ['53', '54', '55', '56'],
    #     "setting1": {
    #         "ood_type": ['Melanocyte', 'Schwann Cell', 'Sebaceous Gland'],
    #         "split_type": "batch", # "cell_type" or "batch"
    #         "train_batch": ['53', '54'],
    #         "test_batch": ['55', '56']}
    # },
    "Human-Kidney-RNA": {
        "modalities": ["RNA", "ATAC"],
        "batches": ['e0def004-9e30-4a3b-9a65-007110f3a1f2', 'f6c0f811-2fb8-4989-b796-37c14b055517', '5028f75a-8c09-4155-a232-ad7dbfa6042e', '8c570254-4bef-48d8-bd79-c812f60835a5', '8213a3f7-2437-4e8a-b836-caec33df901d'],
        "setting3": {
            "info": "",
            "cell_type_key": "cell_type",
            "batch_key": "batch",
            "ood_type": ['HCT116', 'HeLa-S3', 'K562'],  # 'acinar', 'ductal'
            "remove_type": [],
            "train_batch": ['Batch3', 'Batch4'],
            "test_batch": ['Batch1', 'Batch2']
        },
        "setting0": {
            "ood_type": ['CNT', 'ICA', 'TAL', 'DCT', 'PT'],  # 'ICB', 'MES_FIB', 'PT_VCAM1'                     , 'PT_VCAM1'
            "split_type": "batch", # "cell_type" or "batch"
            "train_batch": ['e0def004-9e30-4a3b-9a65-007110f3a1f2', 'f6c0f811-2fb8-4989-b796-37c14b055517', '5028f75a-8c09-4155-a232-ad7dbfa6042e'],
            "test_batch": ['8c570254-4bef-48d8-bd79-c812f60835a5', '8213a3f7-2437-4e8a-b836-caec33df901d']
        },
        "setting2": {
            "ood_type": ['PODO', 'TAL', ],
            "split_type": "batch", # "cell_type" or "batch"
            "train_batch": ['e0def004-9e30-4a3b-9a65-007110f3a1f2', 'f6c0f811-2fb8-4989-b796-37c14b055517', '5028f75a-8c09-4155-a232-ad7dbfa6042e'],
            "test_batch": ['8c570254-4bef-48d8-bd79-c812f60835a5', '8213a3f7-2437-4e8a-b836-caec33df901d']
        },
        "setting1": {
            "cell_type_key": "cell_type",
            "batch_key": "batch",
            "ood_type": ['PODO', 'TAL', 'PEC'],  # PT_VCAM1
            "remove_type": [],
            "split_type": "batch", # "cell_type" or "batch"
            "train_batch": ['e0def004-9e30-4a3b-9a65-007110f3a1f2', 'f6c0f811-2fb8-4989-b796-37c14b055517', '5028f75a-8c09-4155-a232-ad7dbfa6042e'],
            "test_batch": ['8c570254-4bef-48d8-bd79-c812f60835a5', '8213a3f7-2437-4e8a-b836-caec33df901d']
        },
        "setting4": {
            "ood_type": ['ICB', 'ICA', 'PT_VCAM1', 'TAL'],
            "split_type": "batch", # "cell_type" or "batch"
            "train_batch": ['e0def004-9e30-4a3b-9a65-007110f3a1f2', 'f6c0f811-2fb8-4989-b796-37c14b055517', '5028f75a-8c09-4155-a232-ad7dbfa6042e'],
            "test_batch": ['8c570254-4bef-48d8-bd79-c812f60835a5', '8213a3f7-2437-4e8a-b836-caec33df901d']
        },
        "setting5": {
            "ood_type": ['PEC', 'MES_FIB', 'PT_VCAM1', 'TAL', 'PODO'],  # ICB->PEC
            "split_type": "batch", # "cell_type" or "batch"
            "train_batch": ['e0def004-9e30-4a3b-9a65-007110f3a1f2', 'f6c0f811-2fb8-4989-b796-37c14b055517', '5028f75a-8c09-4155-a232-ad7dbfa6042e'],
            "test_batch": ['8c570254-4bef-48d8-bd79-c812f60835a5', '8213a3f7-2437-4e8a-b836-caec33df901d']
        },
        "setting6": {
            "ood_type": ['ICB', 'MES_FIB', 'PT_VCAM1', 'PODO', 'ENDO', 'TAL', 'PT'], # 'TAL' and 'PT' are the top 2 types
            "split_type": "batch", # "cell_type" or "batch"
            "train_batch": ['e0def004-9e30-4a3b-9a65-007110f3a1f2', 'f6c0f811-2fb8-4989-b796-37c14b055517', '5028f75a-8c09-4155-a232-ad7dbfa6042e'],
            "test_batch": ['8c570254-4bef-48d8-bd79-c812f60835a5', '8213a3f7-2437-4e8a-b836-caec33df901d']
        },
        "setting7": {
            "ood_type": ['ICB', 'MES_FIB', 'PT_VCAM1', 'PODO', 'ENDO', 'PC', 'DCT', 'TAL', 'PT'], # 'TAL', 'DCT' and 'PT' are the top 3 types
            "split_type": "batch", # "cell_type" or "batch"
            "train_batch": ['e0def004-9e30-4a3b-9a65-007110f3a1f2', 'f6c0f811-2fb8-4989-b796-37c14b055517', '5028f75a-8c09-4155-a232-ad7dbfa6042e'],
            "test_batch": ['8c570254-4bef-48d8-bd79-c812f60835a5', '8213a3f7-2437-4e8a-b836-caec33df901d']
        },
    },
    "Human-Kidney-ATAC": {
        "modalities": ["RNA", "ATAC"],
        "batches": ['e0def004-9e30-4a3b-9a65-007110f3a1f2', 'f6c0f811-2fb8-4989-b796-37c14b055517', '5028f75a-8c09-4155-a232-ad7dbfa6042e', '8c570254-4bef-48d8-bd79-c812f60835a5', '8213a3f7-2437-4e8a-b836-caec33df901d'],
        "setting1": {
            "cell_type_key": "cell_type",
            "batch_key": "batch",
            "ood_type": ['PODO', 'TAL', 'PEC'],  # PT_VCAM1
            "remove_type": [],
            "split_type": "batch", # "cell_type" or "batch"
            "train_batch": ['e0def004-9e30-4a3b-9a65-007110f3a1f2', 'f6c0f811-2fb8-4989-b796-37c14b055517', '5028f75a-8c09-4155-a232-ad7dbfa6042e'],
            "test_batch": ['8c570254-4bef-48d8-bd79-c812f60835a5', '8213a3f7-2437-4e8a-b836-caec33df901d']
        },
    },
    "Pancreas": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "去掉所有外分泌细胞: acinar(1146个), ductal(1521个)",
            "cell_type_key": "celltype",
            "batch_key": "batch",
            "ood_type": ['beta', 'ductal'],  # 'acinar', 'ductal'
            "remove_type": [],
            "train_batch": ['Baron'],
            "test_batch": ['Segerstolpe']
        },
    },
    "CL_RNA": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_type",
            "batch_key": "batch",
            "ood_type": ['HCT116', 'HeLa-S3', 'K562'],  # 'acinar', 'ductal'
            "remove_type": [],
            "train_batch": ['Batch3', 'Batch4'],
            "test_batch": ['Batch1', 'Batch2']
        },
    },
    "CL_ATAC": {
        "modalities": ["ATAC"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_type",
            "batch_key": "batch",
            "ood_type": ['HCT116', 'HeLa-S3', 'K562'],  # 'acinar', 'ductal'
            "remove_type": [],
            "train_batch": ['Batch3', 'Batch4'],
            "test_batch": ['Batch1', 'Batch2']
        },
    },
    "UP_HK_RNA": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_type",
            "batch_key": "batch",
            "ood_type": ['MES_FIB', 'ICB', 'LEUK'],  #
            "remove_type": [],
            "train_batch": ['8c570254-4bef-48d8-bd79-c812f60835a5', '5028f75a-8c09-4155-a232-ad7dbfa6042e', '8213a3f7-2437-4e8a-b836-caec33df901d'],
            "test_batch": ['e0def004-9e30-4a3b-9a65-007110f3a1f2', 'f6c0f811-2fb8-4989-b796-37c14b055517']
        },
        "setting2": {
            "info": "",
            "cell_type_key": "cell_type",
            "batch_key": "batch",
            "ood_type": ['PEC', 'PODO', 'ENDO'],  #
            "remove_type": [],
            "train_batch": ['8c570254-4bef-48d8-bd79-c812f60835a5', '5028f75a-8c09-4155-a232-ad7dbfa6042e', '8213a3f7-2437-4e8a-b836-caec33df901d'],
            "test_batch": ['e0def004-9e30-4a3b-9a65-007110f3a1f2', 'f6c0f811-2fb8-4989-b796-37c14b055517']
        },
    },
    "UP_HK_ATAC": {
        "modalities": ["ATAC"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_type",
            "batch_key": "batch",
            "ood_type": ['DCT', 'TAL'],  #
            "remove_type": [],
            "train_batch": ['8c570254-4bef-48d8-bd79-c812f60835a5', '5028f75a-8c09-4155-a232-ad7dbfa6042e', '8213a3f7-2437-4e8a-b836-caec33df901d'],
            "test_batch": ['e0def004-9e30-4a3b-9a65-007110f3a1f2', 'f6c0f811-2fb8-4989-b796-37c14b055517']
        },
    },
    "PurifiedPBMCDataset": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_types",
            "batch_key": "batch",
            "ood_type": ['b_cells', 'naive_t', 'naive_cytotoxic'],  # 'acinar', 'ductal'
            "remove_type": [],
            "train_batch": [0, 1, 8, 3, 4, 9, 10, 7],
            "test_batch": [2, 5, 6]
        },
    },
    "TS_Bone_Marrow": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_ontology_class",
            "batch_key": "method",
            "ood_type": ['monocyte', 'cd4-positive, alpha-beta t cell', 'cd8-positive, alpha-beta t cell', 'cd24 neutrophil', 'plasmablast'],  # 'acinar', 'ductal'
            "remove_type": [],
            "train_batch": ['10X'],
            "test_batch": ['smartseq2']
        },
        "setting2": {
            "info": "",
            "cell_type_key": "cell_ontology_class",
            "batch_key": "donor",
            "ood_type": ['plasma cell', 'monocyte'],  # 'acinar', 'ductal'
            "remove_type": ['plasmablast'],
            "train_batch": ['TSP2', 'TSP11', 'TSP13'],
            "test_batch": ['TSP14']
        },
    },
    "TS_Skin": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_ontology_class",
            "batch_key": "method",
            "ood_type": ['macrophage', 'mast cells'],  # 'acinar', 'ductal'
            "remove_type": [],
            "train_batch": ['10X'],
            "test_batch": ['smartseq2']
        },
        "setting2": {
            "info": "",
            "cell_type_key": "cell_ontology_class",
            "batch_key": "donor",
            "ood_type": ['endothelial cell', 'cell of skeletal muscle', 'smooth muscle cell', 'muscle cell', 'stromal cell'],  # 'acinar', 'ductal'
            "remove_type": [],
            "train_batch": ['TSP10'],
            "test_batch": ['TSP14']
        },
    },
    "Lung": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_type",
            "batch_key": "batch",
            "ood_type": ['Macrophage', 'Dendritic cell'],  # 'acinar', 'ductal'
            "remove_type": ['Type 1', 'Basal 1', 'Lymphatic', 'Ionocytes'],
            "train_batch": ['1', '2', '3'],  # , '4', '5'
            "test_batch": ['4', '5', '6']
        },
    },
    "TMS": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "将3/18/21m作为训练集，24/30m作为测试集，ood类为：pancreatic B cell 729个, pancreatic ductal cell 343个",
            "cell_type_key": "cell_ontology_class",
            "batch_key": "age",
            "ood_type": ['pancreatic B cell', 'pancreatic ductal cell'],
            "remove_type": [],
            "train_batch": ['3m', '18m', '21m'],  # '3m',
            "test_batch": ['24m', '30m'],  # , '30m'
        },
        "setting2": {
            "info": "将3/18/21m作为训练集，24/30m作为测试集，ood类为：pancreatic B cell 729个, pancreatic ductal cell 343个",
            "cell_type_key": "cell_ontology_class",
            "batch_key": "age",
            "ood_type": ['endothelial cell', 'leukocyte', 'pancreatic A cell', 'pancreatic D cell', 'pancreatic PP cell'],
            "remove_type": [],
            "train_batch": ['3m'],  # '3m',
            "test_batch": ['18m'],  # , '30m'
        },
    },
    "Mus_Kidney": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_type1",
            "batch_key": "age",
            "ood_type": ['Endothelium', 'Immune', 'Mesangium', 'Tubules'],
            "remove_type": ['B lymph', 'CD-IC', 'CD-PC', 'CD-Trans', 'DCT', 'Distal tubule', 'ED', 'Endo', 'Endothelium', 'Fib', 'Immune', 'LOH', 'Macro', 'NK', 'Neutro', 'Novel1', 'Novel2', 'Podo', 'Podocytes', 'T lymph', 'Tubules', 'endothelial cell', 'epithelial cell of proximal tubule', 'kidney capillary endothelial cell', 'kidney cell', 'kidney collecting duct epithelial cell', 'kidney loop of Henle ascending limb epithelial cell', 'kidney proximal straight tubule epithelial cell', 'macrophage', 'mesangial cell', 'natural killer cell',],
            "train_batch": ['post natal day 1'],  # '3m',
            "test_batch": ['8 weeks'],  # , '30m'
        },
    },
    "Giraddi_10x": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_type1",
            "batch_key": "lifestage",
            "ood_type": ['Luminal trajectory'],
            "remove_type": [],
            "train_batch": ['E16', 'E18'],  # '3m',
            "test_batch": ['Adu1', 'Adu2'],  # , '30m'
        },
    },
    "Han": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_type1",
            "batch_key": "Stages",
            "ood_type": ['Neural_Crest'],
            "remove_type": [],
            "train_batch": ['E8.5'],  # '3m',
            "test_batch": ['E9.0'],  # , '30m'
        },
    },
    "Hochane": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_type1",
            "batch_key": "age",
            "ood_type": ['CnT', 'Leu', 'Mes', 'NPCb', 'Prolif'],
            "remove_type": [],
            "train_batch": ['9 weeks'],  # '3m',
            "test_batch": ['16 weeks'],  # , '30m'
        },
    },
    "Celegans": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell.type",
            "batch_key": "time.point",  # time.point  embryo.time.bin
            "ood_type": ['Excretory_gland', 'XXX'],
            "remove_type": [],  # 'ABarpaaa_lineage', 'Excretory_cell_parent', 'Excretory_duct_and_pore', 'G2_and_W_blasts', 'GLR', 'Parent_of_exc_duct_pore_DB_1_3', 'Parent_of_exc_gland_AVK', 'Parent_of_hyp1V_and_ant_arc_V', 'Pharyngeal_marginal_cell', 'Pharyngeal_neuron', 'T', 'hmc_and_homolog', 'hmc_homolog', 'hyp1V_and_ant_arc_V'
            "train_batch": ['300_minutes'],  # 300_minutes  580-650
            "test_batch": ['500_minutes'],  # 500_minutes  > 650
        },
    },
    "PBMC_Time": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "celltype.l1",
            "batch_key": "time",
            "ood_type": ['other', 'NK', 'other T'],  # , 'other T'
            "remove_type": [],
            "train_batch": ['0', '3'],  #
            "test_batch": ['7'],  #
        },
        "setting2": {
            "info": "",
            "cell_type_key": "celltype.l2",
            "batch_key": "time",
            "ood_type": ['CD4 TCM', 'NK', 'CD8 TEM'],  # , 'other T'
            "remove_type": [],
            "train_batch": ['0', '3'],  #
            "test_batch": ['7'],  #
        },
    },
    "PBMC_Sub": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "sub_cell_type",
            "batch_key": "time",
            "ood_type": ['B naive', 'B memory', 'B intermediate'],  # , 'Plasmablast'
            "remove_type": [],
            "train_batch": ['0', '3'],  # '3m',
            "test_batch": ['7'],  # , '30m'
        },
        "setting2": {
            "info": "",
            "cell_type_key": "sub_cell_type_CD4",
            "batch_key": "time",
            "ood_type": ['CD4 TCM', 'CD4 TEM', 'CD4 CTL', 'CD4 Naive', 'Treg'],  # , 'Plasmablast'
            "remove_type": [],
            "train_batch": ['0', '3'],  # '3m',
            "test_batch": ['7'],  # , '30m'
        },
        "setting3": {
            "info": "",
            "cell_type_key": "sub_cell_type_NK",
            "batch_key": "time",
            "ood_type": ['NK_CD56bright', 'NK'],  # , 'Plasmablast'
            "remove_type": [],
            "train_batch": ['0', '3'],  # '3m',
            "test_batch": ['7'],  # , '30m'
        },
        "setting4": {
            "info": "",
            "cell_type_key": "sub_cell_type_CD8",
            "batch_key": "time",
            "ood_type": ['CD8 TCM', 'CD8 Naive', 'CD8 TEM'],  # , 'Plasmablast'
            "remove_type": [],
            "train_batch": ['0', '3'],  # '3m',
            "test_batch": ['7'],  # , '30m'
        },
        "setting5": {
            "info": "",
            "cell_type_key": "sub_cell_type_B_l3",
            "batch_key": "time",
            "ood_type": ['B intermediate kappa', 'B memory lambda', 'B intermediate lambda', 'B memory kappa', 'B naive lambda', 'B naive kappa'],  # , 'Plasmablast'
            "remove_type": [],
            "train_batch": ['0', '3'],  # '3m',
            "test_batch": ['7'],  # , '30m'
        },
    },
    "BGI": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "celltype_coarse",
            "batch_key": "time_point",
            "ood_type": ['Ectoderm', 'Mesoderm'],
            "remove_type": [],
            "train_batch": ['E7.5', 'E7.75'],  # 'E7.5',
            "test_batch": ['E8.0'],  #
        },
    },
    "Dorsal_Midbrain": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "annotation",
            "batch_key": "Time point",
            "ood_type": ['RGC', 'Glu Neu'],
            "remove_type": [],
            "train_batch": ['E12.5', 'E14.5'],  # , 'E14.5'
            "test_batch": ['E16.5'],
        },
    },
    "human_cortex_EaFet_Adol": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "Cell_type",
            "batch_key": "Donor_ID",
            "ood_type": ['Oligodendrocytes', 'VSMC'],
            "remove_type": [],
            "train_batch": ['EaFet1', 'EaFet2'],  # , 'E14.5'
            "test_batch": ['Adol1', 'Adol2'],
        },
    },
    "Zeisel_p20_p23": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_type1",
            "batch_key": "age",
            "ood_type": ['Dentate gyrus granule neurons', 'Dentate gyrus radial glia-like cells', 'Telencephalon inhibitory interneurons', 'Telencephalon projecting excitatory neurons', 'Non-glutamatergic neuroblasts'],
            "remove_type": [],  # 'Enteric glia', 'Enteric neurons', 'Olfactory ensheathing cells', 'Olfactory inhibitory neurons', 'Sympathetic cholinergic neurons', 'Sympathetic noradrenergic neurons', 'Telencephalon projecting inhibitory neurons'
            "train_batch": ['p20'],  #
            "test_batch": ['p23'],
        },
    },
    "FL_rLN": {
        "modalities": ["ADT"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_type",
            "batch_key": "Entity",
            "ood_type": ['Tfh'],
            "remove_type": [],  # 'Enteric glia', 'Enteric neurons', 'Olfactory ensheathing cells', 'Olfactory inhibitory neurons', 'Sympathetic cholinergic neurons', 'Sympathetic noradrenergic neurons', 'Telencephalon projecting inhibitory neurons'
            "train_batch": ['rLN'],  #
            "test_batch": ['FL'],
        },
    },
    "Zeisel_p20_p29": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_type1",
            "batch_key": "age",
            "ood_type": ['Olfactory ensheathing cells', 'Olfactory inhibitory neurons', 'Dentate gyrus radial glia-like cells', 'Telencephalon inhibitory interneurons', 'Telencephalon projecting excitatory neurons', 'Telencephalon projecting inhibitory neurons', 'Non-glutamatergic neuroblasts'],
            "remove_type": [],  # 'Enteric glia', 'Enteric neurons', 'Olfactory ensheathing cells', 'Olfactory inhibitory neurons', 'Sympathetic cholinergic neurons', 'Sympathetic noradrenergic neurons', 'Telencephalon projecting inhibitory neurons'
            "train_batch": ['p20'],  #
            "test_batch": ['p29'],
        },
        "setting2": {
            "info": "",
            "cell_type_key": "cell_type1",
            "batch_key": "age",
            "ood_type": ['Oligodendrocytes', 'Olfactory ensheathing cells', 'Non-glutamatergic neuroblasts'],
            "remove_type": ['Cerebellum neurons', 'Cholinergic and monoaminergic neurons', 'Choroid epithelial cells', 'Dentate gyrus radial glia-like cells', 'Di- and mesencephalon excitatory neurons', 'Di- and mesencephalon inhibitory neurons', 'Glutamatergic neuroblasts', 'Hindbrain neurons', 'Microglia', 'Olfactory inhibitory neurons', 'Peptidergic neurons', 'Subventricular zone radial glia-like cells', 'Peripheral sensory neurofilament neurons', 'Peripheral sensory non-peptidergic neurons', 'Peripheral sensory peptidergic neurons', 'Perivascular macrophages', 'Satellite glia', 'Schwann cells', 'Spinal cord excitatory neurons', 'Spinal cord inhibitory neurons', 'Subcommissural organ hypendymal cells', 'Telencephalon inhibitory interneurons', 'Telencephalon projecting excitatory neurons', 'Telencephalon projecting inhibitory neurons', 'Pericytes', 'Vascular endothelial cells', 'Vascular smooth muscle cells' ],  # 'Enteric glia', 'Enteric neurons', 'Olfactory ensheathing cells', 'Olfactory inhibitory neurons', 'Sympathetic cholinergic neurons', 'Sympathetic noradrenergic neurons', 'Telencephalon projecting inhibitory neurons'
            "train_batch": ['p20'],  #
            "test_batch": ['p29'],
        },
    },
    "PancreasHumanMouse": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "将老鼠作为训练集，人类作为测试集，ood类为：acinar（腺泡细胞 / 腺泡上皮细胞）958个, mast（肥大细胞）25个，epsilon（ε 细胞 / ε型内分泌细胞）18个",
            "cell_type_key": "cell_type",
            "batch_key": "species",
            "ood_type": ['acinar', 'mast', 'epsilon'],  # , 'mast', 'epsilon'
            "remove_type": ['mast', 'epsilon'],
            "train_batch": ['Mouse'],
            "test_batch": ['Human'],
        },
    },
    "Baron_Human_Mouse_Normalized": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "将老鼠作为训练集，人类作为测试集，ood类为：acinar（腺泡细胞 / 腺泡上皮细胞）958个, mast（肥大细胞）25个，epsilon（ε 细胞 / ε型内分泌细胞）18个",
            "cell_type_key": "cell_type",
            "batch_key": "species",
            "ood_type": ['acinar'],  # , 'mast', 'epsilon'
            "remove_type": ['mast', 'epsilon'],  # 'mast', 'epsilon'
            "train_batch": ['mouse'],
            "test_batch": ['human'],
        },
    },
    "Retina_Human_Mouse": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_type1",
            "batch_key": "species",
            "ood_type": ['Astrocytes', 'Endothelium', 'Ganglion', 'Horizontal', 'Microglia', 'Pericytes'],  #
            "remove_type": [],  #
            "train_batch": ['mouse'],
            "test_batch": ['human'],
        },
    },
    "Retina_Human_Mouse_v2": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_type1",
            "batch_key": "species",
            "ood_type": ['Endothelium'],  #
            "remove_type": [],  #
            "train_batch": ['mouse'],
            "test_batch": ['human'],
        },
    },
    "Retina_Human_Mouse_v3": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_ontology_class",
            "batch_key": "species",
            "ood_type": ['microglial cell', 'retinal ganglion cell'],  #
            "remove_type": [],  #
            "train_batch": ['mouse'],
            "test_batch": ['human'],
        },
    },
    "Thymus_Cao_ALIGNED": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_ontology_class",
            "batch_key": "species",
            "ood_type": ['blood vessel endothelial cell', 'epithelial cell of thymus', 'stromal cell', 'thymocyte'],  #
            "remove_type": ['immature T cell'],  #
            "train_batch": ['mouse'],
            "test_batch": ['human'],
        },
        "setting2": {
            "info": "",
            "cell_type_key": "cell_ontology_class",
            "batch_key": "dataset",
            "ood_type": ['immature T cell'],  #
            "remove_type": [],  #
            "train_batch": ['Mus_musculus'],
            "test_batch": ['Smartseq2'],
        },
    },
    "ImmuneHumanMouse": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "将人作为训练集，鼠作为测试集，ood类为：Neutrophils（中性粒细胞）21928个, Lymphocyte progenitors（淋巴细胞祖细胞7048个",
            "cell_type_key": "final_annotation",
            "batch_key": "species",
            "ood_type": ['Basophils', 'Lymphocyte progenitors', 'Neutrophils'],
            "remove_type": [],
            "train_batch": ['Human'],
            "test_batch": ['Mouse'],
        },
        "setting2": {
            "info": "将人作为训练集，鼠作为测试集，ood类为：Neutrophils（中性粒细胞）21928个, Lymphocyte progenitors（淋巴细胞祖细胞7048个",
            "cell_type_key": "cell_type",
            "batch_key": "species",
            "ood_type": ['B cells', 'Monocytes cells', 'T/NK cells'],
            "remove_type": [],
            "train_batch": ['Mouse'],
            "test_batch": ['Human'],
        },
        "setting3": {
            "info": "将人作为训练集，鼠作为测试集，ood类为：Neutrophils（中性粒细胞）21928个, Lymphocyte progenitors（淋巴细胞祖细胞7048个",
            "cell_type_key": "final_annotation",
            "batch_key": "species",
            "ood_type": ['CD4+ T cells', 'Plasma cells', 'Plasmacytoid dendritic cells'],
            "remove_type": [],
            "train_batch": ['Mouse'],
            "test_batch": ['Human'],
        },
    },
    "Immune_PBMC_Human_Mouse": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "final_annotation",
            "batch_key": "species",
            "ood_type": ['CD4+ T cells', 'NK cells', 'Monocytes', 'CD8+ T cells', 'Plasmacytoid dendritic cells'],
            "remove_type": [],
            "train_batch": ['Mouse'],
            "test_batch": ['Human'],
        },
    },
    "Liver_Human_Monkey": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_ontology_class",
            "batch_key": "species",
            "ood_type": ['Macrophages'],
            "remove_type": [],
            "train_batch": ['monkey'],
            "test_batch": ['human'],
        },
        "setting2": {
            "info": "只有一类重合",
            "cell_type_key": "cell_ontology_class",
            "batch_key": "species",
            "ood_type": ['Macrophages', 'Endothelial', 'fibroblasts', 'Cholangio'],  # , 'Endothelial' , 'Cholangio'  'Macrophages',  , 'fibroblasts'
            "remove_type": [],
            "train_batch": ['monkey'],
            "test_batch": ['human'],
        },
        "setting3": {
            "info": "没有类别重合",
            "cell_type_key": "cell_ontology_class",
            "batch_key": "species",
            "ood_type": ['Macrophages', 'Cholangio', 'Endothelial', 'fibroblasts', 'Hepatocytes'],
            "remove_type": [],
            "train_batch": ['monkey'],
            "test_batch": ['human'],
        },
        "setting4": {
            "info": "",
            "cell_type_key": "cell_ontology_class",
            "batch_key": "species",
            "ood_type": ['B Plasma cells', 'Kupffer cells', 'NK NKT and T cells'],
            "remove_type": [],
            "train_batch": ['human'],
            "test_batch": ['monkey'],
        },
    },
    "Liver_Human_Zebrafish": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_ontology_class",
            "batch_key": "species",
            "ood_type": ['Cholangio', 'fibroblasts'],
            "remove_type": [],
            "train_batch": ['zebrafish'],
            "test_batch": ['human'],
        },
    },
    "Liver_Human_Mouse": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_ontology_class",
            "batch_key": "species",
            "ood_type": ['Macrophages'],  # , 'fibroblasts'
            "remove_type": [],
            "train_batch": ['mouse'],
            "test_batch": ['human'],
        },
        "setting2": {
            "info": "",
            "cell_type_key": "cell_ontology_class",
            "batch_key": "species",
            "ood_type": ['B Plasma cells', 'HsPCs', 'Kupffer cells'],
            "remove_type": [],
            "train_batch": ['human'],
            "test_batch": ['mouse'],
        },
    },
    "Testis_Human_Platypus": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_type",  # cell_ontology_class
            "batch_key": "species",
            "ood_type": ['Differentiated_Spermatogonia', 'Early_Round_spermatids', 'Late_Round_spermatids', 'Leptotene_Spermatocytes', 'Pachytene_Spermatocytes', 'Sertoli', 'Undifferentiated_Spermatogonia', 'Zygotene_Spermatocytes'],
            "remove_type": [],
            "train_batch": ['platypus'],
            "test_batch": ['human'],
        },
    },
    "Brain_Human_Turtle": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_ontology_class",
            "batch_key": "species",
            "ood_type": ['Purkinje cell', 'cerebellar granule cell', 'endothelial cell'],
            "remove_type": [],
            "train_batch": ['turtle'],
            "test_batch": ['human'],
        },
    },
    "Brain_Human_Mouse": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_ontology_class",
            "batch_key": "species",
            "ood_type": ['Purkinje cell', 'cerebellar granule cell', 'brain pericyte'],
            "remove_type": [],
            "train_batch": ['mouse'],
            "test_batch": ['human'],
        },
    },
    "Glia_Human_Rhesus": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_ontology_class",
            "batch_key": "species",
            "ood_type": ['Micro P2RY12 CCL3', 'Micro P2RY12 GLDN'],
            "remove_type": [],
            "train_batch": ['rhesus'],
            "test_batch": ['human'],
        },
    },
    "Colorectal_Cancer": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "standard_true_celltype",
            "batch_key": "split",
            "ood_type": ['Cancer_cells'],  # , 'other T'
            "remove_type": [],
            "train_batch": ['train'],  #
            "test_batch": ['test'],  #
        },
        "setting2": {
            "info": "",
            "cell_type_key": "cell_type_stage",
            "batch_key": "split",
            "ood_type": ['I', 'II', 'III', 'IV'],  # , 'other T'
            "remove_type": [],
            "train_batch": ['train'],  #
            "test_batch": ['test'],  #
        },
    },
    "Breast_Cancer": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "standard_true_celltype",
            "batch_key": "split",
            "ood_type": ['Cancer_cells'],  # , 'other T'
            "remove_type": [],
            "train_batch": ['test'],  #
            "test_batch": ['train'],  #
        },
    },
    "RCC": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "standard_true_celltype",
            "batch_key": "split",
            "ood_type": ['Cancer_cells'],  # , 'other T'
            "remove_type": [],
            "train_batch": ['train'],  #
            "test_batch": ['test'],  #
        },
    },
    "HGSOC": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "standard_true_celltype",
            "batch_key": "split",
            "ood_type": ['Cancer_cells'],  # , 'other T'
            "remove_type": [],
            "train_batch": ['train'],  #
            "test_batch": ['test'],  #
        },
    },
    "BMMC_Kidney": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_type",
            "batch_key": "Tissue",
            "ood_type": ['PT', 'TAL', 'DCT', 'CNT', 'ICA', 'PC', 'ENDO', 'PEC', 'PODO', 'PT_VCAM1', 'MES_FIB', 'ICB', 'LEUK'],  # , 'other T'
            "remove_type": [],
            "train_batch": ['BMMC'],  #
            "test_batch": ['UP_HK'],  #
        },
    },
    "mice_zebrafish": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_type",
            "batch_key": "tissue",
            "ood_type": ['exocrine cell', 'type B pancreatic cell', 'pancreatic A cell', 'pancreatic D cell', 'hepatic stellate cell', 'pancreatic PP cell', 'Schwann cell'],  # , 'other T'
            "remove_type": ['Kupffer cell', 'hepatocyte', 'neutrophil', 'mature NK T cell', 'kidney tubule cell', 'unknown'],
            "train_batch": ['heart'],  #
            "test_batch": ['pancreas'],  #
        },
    },
    "mouse_heart_pancreas": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_type_disease",
            "batch_key": "tissue",
            "ood_type": ['exocrine cell', 'exocrine cell (myocardial infarction)', 'type B pancreatic cell', 'type B pancreatic cell (myocardial infarction)', 'pancreatic A cell', 'pancreatic A cell (myocardial infarction)', 'pancreatic D cell', 'pancreatic D cell (myocardial infarction)', 'hepatic stellate cell', 'hepatic stellate cell (myocardial infarction)', 'pancreatic PP cell', 'pancreatic PP cell (myocardial infarction)', 'Schwann cell', 'Schwann cell (myocardial infarction)'],  # , 'other T'
            "remove_type": ['Kupffer cell', 'hepatocyte', 'neutrophil', 'mature NK T cell', 'kidney tubule cell', 'unknown'],
            "train_batch": ['heart'],  #
            "test_batch": ['pancreas'],  #
        },
    },
    "Mouse_Pancreas_Subtissue": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_ontology_class",
            "batch_key": "subtissue",
            "ood_type": ['pancreatic A cell', 'pancreatic B cell', 'pancreatic D cell', 'pancreatic PP cell', '', '', ''],
            "remove_type": [],
            "train_batch": ['Exocrine'],  #
            "test_batch": ['Endocrine'],  #
        },
    },
    "TMS_Fat_Limb": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_ontology_class",
            "batch_key": "tissue",
            "ood_type": ['NK cell', 'epithelial cell'],  # , 'other T'
            "remove_type": [],
            "train_batch": ['Limb'],  #
            "test_batch": ['Fat'],  #
        },
        "setting2": {
            "info": "",
            "cell_type_key": "cell_ontology_class",
            "batch_key": "tissue",
            "ood_type": ['skeletal muscle satellite cell'],  # , 'other T'
            "remove_type": [],
            "train_batch": ['Fat'],  #
            "test_batch": ['Limb'],  #
        },
    },
    "TMS_marrow_spleen": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_ontology_class",
            "batch_key": "tissue",
            "ood_type": ['granulocyte_lineage', 'hematopoietic_stem_progenitor', 'myeloid_progenitor'],  # , 'other T'
            "remove_type": [],
            "train_batch": ['Spleen'],  #
            "test_batch": ['Marrow'],  #
        },
    },
    "ALIGNED_Mus_musculus_Mammary_Gland": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_ontology_class",
            "batch_key": "dataset_name",
            "ood_type": ['B cell', 'T cell', 'macrophage', 'endothelial cell', 'stromal cell'],  # , 'other T'
            "remove_type": [],  # 'mammary alveolar cell', 'myoepithelial cell of mammary gland'
            "train_batch": ['Bach'],  #
            "test_batch": ['Quake_10x_Mammary_Gland'],  #
        },
    },
    "TS_Prostate": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_ontology_class",
            "batch_key": "donor",
            "ood_type": ['b cell', 'cd4-positive, alpha-beta t cell', 'cd8-positive, alpha-beta t cell', 'double-positive, alpha-beta thymocyte', 'fibroblast', 'ionocyte', 'mast cell', 'mucus secreting cell', 'neutrophil', 'serous cell of epithelium of trachea', 'smooth muscle cell', 't cell', 'tracheal goblet cell'],  # , 'other T'
            "remove_type": [],  # 'mammary alveolar cell', 'myoepithelial cell of mammary gland'
            "train_batch": ['TSP2'],  #
            "test_batch": ['TSP6'],  #
        },
    },
    "TMS_Thymus": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_ontology_class",
            "batch_key": "sex",
            "ood_type": ['double negative T cell', 'immature T cell'],  # , 'other T'
            "remove_type": [],  # 'mammary alveolar cell', 'myoepithelial cell of mammary gland'
            "train_batch": ['male'],  #
            "test_batch": ['female'],  #
        },
    },
    "Mouse-MOp-RNA": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_type",
            "batch_key": "Seq_batch",
            "ood_type": ['NP', 'Sst', 'Vip'],  # , 'other T'
            "remove_type": [],  # 'mammary alveolar cell', 'myoepithelial cell of mammary gland'
            "train_batch": ['RTX-841', 'RTX-842'],  #
            "test_batch": ['RTX-850'],  #
        },
    },
    "TS_Eye": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_ontology_class",
            "batch_key": "method",
            "ood_type": ['adipocyte', 'corneal epithelial cell', 'ciliary body', 'epithelial cell of lacrimal sac', 'limbal stem cell', 'ocular surface cell', 'retina horizontal cell', 'retinal ganglion cell'],  # , 'other T'
            "remove_type": [],  # 'mammary alveolar cell', 'myoepithelial cell of mammary gland'
            "train_batch": ['smartseq2'],  # 'CEMBA171206_3C', 'CEMBA171207_3C', , 'CEMBA171212_4B', 'CEMBA171213_4B', 'CEMBA180104_4B'
            "test_batch": ['10X'],  # 'CEMBA180409_2C',
        },
    },
    "TS_Pancreas": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_ontology_class",
            "batch_key": "donor",
            "ood_type": ['pancreatic acinar cell', 'pancreatic ductal cell', 'b cell', 'pancreatic alpha cell', 'pancreatic delta cell', 'pancreatic pp cell'],  # , 'other T'
            "remove_type": [],  # 'mammary alveolar cell', 'myoepithelial cell of mammary gland'
            "train_batch": ['TSP9'],  # 'CEMBA171206_3C', 'CEMBA171207_3C', , 'CEMBA171212_4B', 'CEMBA171213_4B', 'CEMBA180104_4B'
            "test_batch": ['TSP1'],  # 'CEMBA180409_2C',
        },
    },
    "TMS-facs-SCAT": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_ontology_class",
            "batch_key": "age",
            "ood_type": ['mesenchymal stem cell of adipose', 'myeloid cell', 'epithelial cell'],  # , 'other T'
            "remove_type": [],  # 'mammary alveolar cell', 'myoepithelial cell of mammary gland'
            "train_batch": ['3m', '18m'],  #
            "test_batch": ['24m'],  #
        },
    },
    "TMS-facs-Trachea": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_ontology_class",
            "batch_key": "age",
            "ood_type": ['macrophage', 'fibroblast', 'T cell'],  # , 'other T'
            "remove_type": [],  # 'mammary alveolar cell', 'myoepithelial cell of mammary gland'
            "train_batch": ['3m', '18m'],  #
            "test_batch": ['24m'],  #
        },
    },
    "TMS-facs-MAT": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_ontology_class",
            "batch_key": "age",
            "ood_type": ['mesenchymal stem cell of adipose', 'B cell', 'myeloid cell'],  # , 'other T'
            "remove_type": [],  # 'mammary alveolar cell', 'myoepithelial cell of mammary gland'
            "train_batch": ['3m', '18m'],  #
            "test_batch": ['24m'],  #
        },
    },
    "TMS-droplet-Pancreas": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_ontology_class",
            "batch_key": "age",
            "ood_type": ['pancreatic acinar cell', 'pancreatic B cell', 'pancreatic ductal cel'],  # , 'other T'
            "remove_type": [],  # 'mammary alveolar cell', 'myoepithelial cell of mammary gland'
            "train_batch": ['3m', '18m', '24m'],  # 'CEMBA171206_3C', 'CEMBA171207_3C', , 'CEMBA171212_4B', 'CEMBA171213_4B', 'CEMBA180104_4B'
            "test_batch": ['30m'],  # 'CEMBA180409_2C',
        },
    },
    "TMS-facs-Limb_Muscle": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_ontology_class",
            "batch_key": "age",
            "ood_type": ['endothelial cell', 'mesenchymal stem cell'],  # , 'other T'
            "remove_type": [],  # 'mammary alveolar cell', 'myoepithelial cell of mammary gland'
            "train_batch": ['3m', '18m'],  # 'CEMBA171206_3C', 'CEMBA171207_3C', , 'CEMBA171212_4B', 'CEMBA171213_4B', 'CEMBA180104_4B'
            "test_batch": ['24m'],  # 'CEMBA180409_2C',
        },
    },
    "TMS-facs-Lung": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_ontology_class",
            "batch_key": "age",
            "ood_type": ['myeloid dendritic cell', 'fibroblast of lung', 'B cell', 'mature NK T cell'],  # , 'other T'
            "remove_type": [],  # 'mammary alveolar cell', 'myoepithelial cell of mammary gland'
            "train_batch": ['3m', '18m'],  # 'CEMBA171206_3C', 'CEMBA171207_3C', , 'CEMBA171212_4B', 'CEMBA171213_4B', 'CEMBA180104_4B'
            "test_batch": ['24m'],  # 'CEMBA180409_2C',
        },
    },
    "TMS-facs-Mammary_Gland": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_ontology_class",
            "batch_key": "age",
            "ood_type": ['basal cell', 'luminal epithelial cell of mammary gland'],  # , 'other T'
            "remove_type": [],  # 'mammary alveolar cell', 'myoepithelial cell of mammary gland'
            "train_batch": ['3m', '18m'],  # 'CEMBA171206_3C', 'CEMBA171207_3C', , 'CEMBA171212_4B', 'CEMBA171213_4B', 'CEMBA180104_4B'
            "test_batch": ['24m'],  # 'CEMBA180409_2C',
        },
    },
    "TMS-facs-Spleen": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_ontology_class",
            "batch_key": "sex",
            "ood_type": ['CD4-positive, alpha-beta T cell', 'CD8-positive, alpha-beta T cell'],  # , 'other T'
            "remove_type": [],  # 'mammary alveolar cell', 'myoepithelial cell of mammary gland'
            "train_batch": ['male'],  # 'CEMBA171206_3C', 'CEMBA171207_3C', , 'CEMBA171212_4B', 'CEMBA171213_4B', 'CEMBA180104_4B'
            "test_batch": ['female'],  # 'CEMBA180409_2C',
        },
    },
    "qian": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "standard_true_celltype",
            "batch_key": "split",
            "ood_type": ['Cancer_cells'],  #
            "remove_type": [],  #
            "train_batch": ['train'],  #
            "test_batch": ['test'],  #
        },
    },
    "peng": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "standard_true_celltype",
            "batch_key": "split",
            "ood_type": ['Cancer_cells'],  #
            "remove_type": [],  #
            "train_batch": ['train'],  #
            "test_batch": ['test'],  #
        },
    },
    "nath": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "standard_true_celltype",
            "batch_key": "split",
            "ood_type": ['Cancer_cells'],  #
            "remove_type": [],  #
            "train_batch": ['train'],  #
            "test_batch": ['test'],  #
        },
    },
    "bi": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "standard_true_celltype",
            "batch_key": "split",
            "ood_type": ['Cancer_cells'],  #
            "remove_type": [],  #
            "train_batch": ['train'],  #
            "test_batch": ['test'],  #
        },
    },
    "bassez": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "standard_true_celltype",
            "batch_key": "split",
            "ood_type": ['Cancer_cells'],  #
            "remove_type": [],  #
            "train_batch": ['train'],  #
            "test_batch": ['test'],  #
        },
    },
    "Kidney_HeartA": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "CellType",
            "batch_key": "tissue",
            "ood_type": ['Cardiomyocytes', 'Immature B cells', 'Microglia', 'Oligodendrocytes'],  #
            "remove_type": [],  #
            "train_batch": ['Kidney'],  #
            "test_batch": ['Heart'],  #
        },
    },
    "Kidney_Liver": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "CellType",
            "batch_key": "tissue",
            "ood_type": ['Alveolar macrophages', 'Collecting duct', 'DCT/CD', 'Distal convoluted tubule', 'Hematopoietic progenitors', 'Loop of henle', 'Proximal tubule', 'Proximal tubule S3', 'Sperm', 'Type II pneumocytes'],  #
            "remove_type": [],  #
            "train_batch": ['Liver'],  #
            "test_batch": ['Kidney'],  #
        },
    },
    "BoneMarrowA_BoneMarrowB": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "CellType",
            "batch_key": "Batch",
            "ood_type": ['Erythroblasts', 'Monocytes', 'Immature B cells'],  #
            "remove_type": [],  #
            "train_batch": ['BoneMarrowB'],  #
            "test_batch": ['BoneMarrowA'],  #
        },
    },
    "MosA1_MosA2": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "CellType",
            "batch_key": "Batch",
            "ood_type": ['Astrocytes', 'Endothelial Cells', 'Inhibitory Neurons'],  #
            "remove_type": [],  #
            "train_batch": ['MosA1'],  #
            "test_batch": ['MosA2'],  #
        },
    },
    "LargeIntestineA_LargeIntestineB": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "CellType",
            "batch_key": "Batch",
            "ood_type": ['B cells', 'T cells', 'Unknown', 'Microglia', 'Sperm', 'Enterocytes'],  #
            "remove_type": [],  #
            "train_batch": ['LargeIntestineB'],  #
            "test_batch": ['LargeIntestineA'],  #
        },
    },
    "LungA_LungB": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "CellType",
            "batch_key": "Batch",
            "ood_type": ['Endothelial II cells', 'Endothelial I cells', 'Type I pneumocytes', 'Type II pneumocytes', 'T cells', 'Unknown', 'Cardiomyocytes', 'Sperm', 'B cells'],  #
            "remove_type": [],  #
            "train_batch": ['LungB'],  #
            "test_batch": ['LungA'],  #
        },
    },
    "MosM1_MosM2": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "CellType",
            "batch_key": "Batch",
            "ood_type": ['Astrocytes', 'Inhibitory Neurons', 'Microglia'],  #
            "remove_type": [],  #
            "train_batch": ['MosM1'],  #
            "test_batch": ['MosM2'],  #
        },
    },
    "MosP1_MosP2": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "CellType",
            "batch_key": "Batch",
            "ood_type": ['Astrocytes', 'Inhibitory Neurons', 'Microglia'],  #
            "remove_type": [],  #
            "train_batch": ['MosP2'],  #
            "test_batch": ['MosP1'],  #
        },
    },
    "WholeBrainA_WholeBrainB": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "CellType",
            "batch_key": "Batch",
            "ood_type": ['Cerebellar granule cells', 'Inhibitory neurons', 'Astrocytes', 'B cells', 'Regulatory T cells', 'Sperm'],  #
            "remove_type": [],  #
            "train_batch": ['WholeBrainB'],  #
            "test_batch": ['WholeBrainA'],  #
        },
    },
    "TMS-facs-GAT": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_ontology_class",
            "batch_key": "age",
            "ood_type": ['mesenchymal stem cell of adipose', 'B cell', 'myeloid cell'],  #
            "remove_type": [],  #
            "train_batch": ['3m', '18m'],  #
            "test_batch": ['24m'],  #
        },
    },
    "TMS-facs-Liver": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_ontology_class",
            "batch_key": "age",
            "ood_type": ['hepatocyte', 'Kupffer cell', 'myeloid leukocyte'],  #
            "remove_type": [],  #
            "train_batch": ['3m', '18m'],  #
            "test_batch": ['24m'],  #
        },
    },
    "TMS-droplet-Large_Intestine": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_ontology_class",
            "batch_key": "subtissue",
            "ood_type": ['enterocyte of epithelium of large intestine', 'large intestine goblet cell'],  #
            "remove_type": [],  #
            "train_batch": ['COLON - "PROM Tm"', 'COLON PROXIMAL'],  #
            "test_batch": ['COLON:P+D'],  #
        },
    },
    "TMS-facs-BAT": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_ontology_class",
            "batch_key": "age",
            "ood_type": ['epithelial cell'],  #
            "remove_type": [],  #
            "train_batch": ['18m'],  #
            "test_batch": ['24m'],  #
        },
    },
    "TMS-facs-Diaphragm": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_ontology_class",
            "batch_key": "age",
            "ood_type": ['skeletal muscle satellite cell', 'endothelial cell'],  #
            "remove_type": [],  #
            "train_batch": ['3m', '18m'],  #
            "test_batch": ['24m'],  #
        },
    },
    "TMS-facs-Kidney": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_ontology_class",
            "batch_key": "sex",
            "ood_type": ['kidney interstitial fibroblast', 'kidney collecting duct principal cell', 'kidney collecting duct epithelial cell', 'epithelial cell of proximal tubule'],  #
            "remove_type": [],  #
            "train_batch": ['female'],  #
            "test_batch": ['male'],  #
        },
    },
    "TMS-facs-Aorta": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_ontology_class",
            "batch_key": "sex",
            "ood_type": ['aortic endothelial cell', 'fibroblast of cardiac tissue'],  #
            "remove_type": [],  #
            "train_batch": ['male'],  #
            "test_batch": ['female'],  #
        },
    },
    "Cao_2020_Eye": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_ontology_class",
            "batch_key": "Batch",
            "ood_type": ['amacrine cell', 'retinal ganglion cell', 'retina horizontal cell'],  #
            "remove_type": ['nan'],  #
            "train_batch": ['2', '4'],  #
            "test_batch": ['10'],  #
        },
    },
    "Hu": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_ontology_class",
            "batch_key": "donor",
            "ood_type": ['Muller cell', 'amacrine cell', 'fibroblast', 'photoreceptor cell', 'retinal bipolar neuron'],  #
            "remove_type": ['nan'],  #
            "train_batch": ['emb3'],  #
            "test_batch": ['emb2'],  #
        },
    },
    "TMS-droplet-Kidney": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_ontology_class",
            "batch_key": "sex",
            "ood_type": ['kidney proximal convoluted tubule epithelial cell', 'fenestrated cell ', 'kidney distal convoluted tubule epithelial cell'],  #
            "remove_type": [],  #
            "train_batch": ['male'],  #
            "test_batch": ['female'],  #
        },
    },
    "TMS-Skin": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_ontology_class",
            "batch_key": "sex",
            "ood_type": ['basal cell of epidermis', 'epidermal cell'],  #
            "remove_type": [],  #
            "train_batch": ['male'],  #
            "test_batch": ['female'],  #
        },
    },
    "cd34_multiome_rna": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "celltype",
            "batch_key": "batch",
            "ood_type": ['HSC', 'Mono', 'Ery'],  #
            "remove_type": [],  #
            "train_batch": ['0'],  #
            "test_batch": ['1'],  #
        },
    },
    "bm_multiome_rna": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "celltype",
            "batch_key": "batch",
            "ood_type": ['Mono', 'Mono2', 'MonoPre'],  #
            "remove_type": [],  #
            "train_batch": ['0'],  #
            "test_batch": ['1'],  #
        },
    },
    "TS_Thymus": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_ontology_class",
            "batch_key": "donor",
            "ood_type": ['b cell', 'cd8-positive, alpha-beta cytotoxic t cell', 'dn1 thymic pro-t cell', 'dn3 thymocyte', 'dn4 thymocyte', 'innate lymphoid cell', 'myeloid dendritic cell', 'naive regulatory t cell', 't follicular helper cell'],  #
            "remove_type": [],  #
            "train_batch": ['TSP2'],  #
            "test_batch": ['TSP14'],  #
        },
    },
    "PBMC_ADT": {
        "modalities": ["ADT"],
        "setting1": {
            "info": "",
            "cell_type_key": "celltype.l2",
            "batch_key": "time",
            "ood_type": ['B naive', 'B memory', 'B intermediate'],  # , 'Plasmablast'
            "remove_type": [],
            "train_batch": ['0', '3'],  # '3m',
            "test_batch": ['7'],  # , '30m'
        },
    },
    "BMMC_RNA": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_type",
            "batch_key": "is_train",
            "ood_type": ['CD14+ Mono', 'Naive CD20+ B IGKC+', 'Naive CD20+ B IGKC-', 'Proerythroblast'],  #
            "remove_type": [],  #
            "train_batch": ['train'],  #
            "test_batch": ['test'],  #
        },
    },
    "BMMC_ADT": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_type",
            "batch_key": "is_train",
            "ood_type": ['CD14+ Mono', 'Naive CD20+ B IGKC+', 'Naive CD20+ B IGKC-', 'Proerythroblast'],  #
            "remove_type": [],  #
            "train_batch": ['train'],  #
            "test_batch": ['test'],  #
        },
    },
    "Ma_RNA": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_type",
            "batch_key": "batch",
            "ood_type": ['Infundibulum', 'Spinous', 'TAC-1', 'TAC-2'],  #
            "remove_type": [],  #
            "train_batch": ['53', '54'],  #
            "test_batch": ['55', '56'],  #
        },
    },
    "Ma_ATAC": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_type",
            "batch_key": "batch",
            "ood_type": ['Infundibulum', 'Spinous', 'TAC-1', 'TAC-2'],  #
            "remove_type": [],  #
            "train_batch": ['53', '54'],  #
            "test_batch": ['55', '56'],  #
        },
    },
    "Cao_2020_Pancreas": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_ontology_class",
            "batch_key": "donor",
            "ood_type": ['blood vessel endothelial cell', 'leukocyte', 'stromal cell of pancreas'],  #
            "remove_type": [],  #
            "train_batch": ['H27870', 'H27876'],  #
            "test_batch": ['H27948'],  #
        },
    },
    "Alemany_Fin": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_ontology_class",
            "batch_key": "donor",
            "ood_type": ['fibroblast', 'leukocyte'],  #
            "remove_type": [],  #
            "train_batch": ['R4'],  #
            "test_batch": ['R5'],  #
        },
    },
    "mice_kidney_heart": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_type",
            "batch_key": "tissue",
            "ood_type": ['fibroblast', 'T cell', 'granulocyte', 'fibrocyte', 'dendritic cell', 'monocyte', 'natural killer cell', 'cardiac muscle cell', 'endocardial cell'],  #
            "remove_type": ['unknown', 'Kupffer cell', 'exocrine cell', 'type B pancreatic cell', 'pancreatic A cell', 'pancreatic D cell', 'hepatocyte', 'hepatic stellate cell', 'pancreatic PP cell', 'Schwann cell'],  #
            "train_batch": ['kidney'],  #
            "test_batch": ['heart'],  #
        },
        "setting2": {
            "info": "",
            "cell_type_key": "cell_type",
            "batch_key": "tissue",
            "ood_type": ['neutrophil', 'mature NK T cell', 'kidney tubule cell'],  #
            "remove_type": ['unknown', 'Kupffer cell', 'exocrine cell', 'type B pancreatic cell', 'pancreatic A cell', 'pancreatic D cell', 'hepatocyte', 'hepatic stellate cell', 'pancreatic PP cell', 'Schwann cell'],  #
            "train_batch": ['heart'],  #
            "test_batch": ['kidney'],  #
        },
    },
    "mice_liver_heart": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_type",
            "batch_key": "tissue",
            "ood_type": ['Kupffer cell', 'hepatocyte', 'hepatic stellate cell'],  #
            "remove_type": ['unknown', 'exocrine cell', 'type B pancreatic cell', 'pancreatic A cell', 'neutrophil', 'mature NK T cell', 'kidney tubule cell',
                            'pancreatic D cell', 'pancreatic PP cell', 'Schwann cell'],  #
            "train_batch": ['heart'],  #
            "test_batch": ['liver'],  #
        },
    },
    "mice_blood_heart": {
        "modalities": ["RNA"],
        "setting1": {
            "info": "",
            "cell_type_key": "cell_type",
            "batch_key": "tissue",
            "ood_type": ['endothelial cell', 'fibrocyte', 'pericyte', 'endocardial cell'],  #
            "remove_type": ['Kupffer cell', 'exocrine cell', 'type B pancreatic cell', 'pancreatic A cell', 'pancreatic D cell', 'hepatic stellate cell', 'neutrophil', 'mature NK T cell', 'pancreatic PP cell', 'Schwann cell', 'kidney tubule cell', 'unknown'],  #
            "train_batch": ['blood'],  #
            "test_batch": ['heart'],  #
        },
    },
    "PBMC_ATAC_ASAP_seq": {
        "modalities": ["ATAC"],
        "setting1": {
            "info": "",
            "cell_type_key": "CellType",
            "batch_key": "split",
            "ood_type": ['Monocytes', 'B', 'NK', 'Effector CD4+ T'],  #
            "remove_type": [],  #
            "train_batch": ['train'],  #
            "test_batch": ['test'],  #
        },
    },
    "MouseKidney_ATAC": {
        "modalities": ["ATAC"],
        "setting1": {
            "info": "",
            "cell_type_key": "CellType",
            "batch_key": "split",
            "ood_type": ['Proximal tubule brush', 'Proximal tubule cell', 'Distal tubule cell', 'Distal convoluted tub', 'Distal collecting duc'],  #
            "remove_type": [],  #
            "train_batch": ['train'],  #
            "test_batch": ['test'],  #
        },
    },
    "Human_Bone_Marrow_CITE_ADT": {
        "modalities": ["ADT"],
        "setting1": {
            "info": "",
            "cell_type_key": "CellType",
            "batch_key": "split",
            "ood_type": ['CD4 Memory', 'CD4 Naive', 'CD14 Mono', 'CD16 Mono'],  #
            "remove_type": [],  #
            "train_batch": ['train'],  #
            "test_batch": ['test'],  #
        },
    },
    "BMMC": {
        "modalities": ["RNA", "ADT"],
        "setting1": ['CD4+ T CD314+ CD45RA+', 'dnT', 'CD8+ T naive CD127+ CD26- CD101-', 'T prog cycling', 'cDC1'],
    },
}