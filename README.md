# scDiscovery

The official implementation for **scDiscovery**.

scDiscovery is a computational framework for single-cell data analysis and discovery. This repository provides code, tutorials, and example workflows for preprocessing single-cell data, learning informative cell representations, and performing downstream discovery tasks such as clustering, annotation, visualization, and prediction.

## Table of Contents

- [Datasets](#datasets)
- [Installation](#installation)
- [Usage](#usage)
- [Tutorial](#tutorial)
- [Output](#output)
- [Citation](#citation)

## Datasets

We provide convenient access to the downloadable datasets through the download links listed in the [`dataset`](./data/dataset.md) text file.

After downloading the datasets, please organize the files as follows:

```text
scDiscovery/
├── data/
│   ├── dataset.txt
|   └── dataset1
│   └── dataset2
|   └── ...
├── output
├── figures
├── scDiscovery.ipynb
├── README.md
└── requirements.txt
```

## Installation

To reproduce scDiscovery, we suggest first creating a conda environment:

```bash
conda create -n scDiscovery python=3.11
conda activate scDiscovery
```

Then clone this repository:

```bash
git clone https://github.com/Hzh-ml/scDiscovery.git
cd scDiscovery
```

Install the required packages:

```bash
pip install -r requirements.txt
```

## Usage

### Data configuration

Before running scDiscovery, prepare the input single-cell dataset in `.h5ad` format and place it in the `data` directory.

Specify the dataset name and its corresponding file in `data_config.py`:

```python
dataset_dict = {
    "example_dataset": "example_dataset_file_name"
}
```

### Training and prediction

After data configuration, use `load_data()` to load, preprocess, and return the processed data.

Run `train_single_omics()` to train scDiscovery and learn cell representations. Then, use `adaptive_decision_boundary_calibration()` to estimate an entropy-based decision boundary for discriminating known and unknown cells.

Finally, call `dynamic_novel_cell_type_expansion()` to predict known cell types and resolve unknown cells into fine-grained novel cell populations. The output `AnnData` object contains the learned scDiscovery embeddings and predicted cell types.

### Evaluation

Use `evaluate_discovery_potential()` and `evaluate_model_on_novel_cell_type()` to evaluate discovery performance.

For more detailed instructions, please refer to and run the tutorial notebook [`scDiscovery.ipynb`](scDiscovery.ipynb).

## Tutorial

### Tutorial 1: Cross-tissue novel cell type discovery (Heart-Pancreas)

1. Install the required environment according to [Installation](#installation).
2. Please download the datasets using the download links provided in the dataset file and place the downloaded files in the `data` directory. If you use your own data, please preprocess it into `.h5ad` format and place the resulting files in the `data` directory.
3. For more detailed instructions, please refer to and run the tutorial notebook [`scDiscovery_cross_tissue.ipynb`](scDiscovery_cross_tissue.ipynb).

### Tutorial 2: Cross-species novel cell type discovery (Liver-Human-Monkey)

1. Install the required environment according to [Installation](#installation).
2. Please download the datasets using the download links provided in the dataset file and place the downloaded files in the `data` directory. If you use your own data, please preprocess it into `.h5ad` format and place the resulting files in the `data` directory.
3. For more detailed instructions, please refer to and run the tutorial notebook [`scDiscovery_cross_species.ipynb`](scDiscovery_cross_species.ipynb).

### Tutorial 3: Cross-developmental-stage novel cell type discovery (Zeisel-2018)

1. Install the required environment according to [Installation](#installation).
2. Please download the datasets using the download links provided in the dataset file and place the downloaded files in the `data` directory. If you use your own data, please preprocess it into `.h5ad` format and place the resulting files in the `data` directory.
3. For more detailed instructions, please refer to and run the tutorial notebook [`scDiscovery_cross_developmental_stages.ipynb`](scDiscovery_cross_developmental_stages.ipynb).

### Tutorial 4: Cancer cell discovery (Peng-PDAC)

1. Install the required environment according to [Installation](#installation).
2. Please download the datasets using the download links provided in the dataset file and place the downloaded files in the `data` directory. If you use your own data, please preprocess it into `.h5ad` format and place the resulting files in the `data` directory.
3. For more detailed instructions, please refer to and run the tutorial notebook [`scDiscovery_cancer_cell_discovery.ipynb`](scDiscovery_cancer_cell_discovery.ipynb).

### Tutorial 5: Novel cell population discovery across distinct cell states (HNSCC-RNA and HNSCC-ADT)

1. Install the required environment according to [Installation](#installation).
2. Please download the datasets using the download links provided in the dataset file and place the downloaded files in the `data` directory. If you use your own data, please preprocess it into `.h5ad` format and place the resulting files in the `data` directory.
3. For more detailed instructions, please refer to and run the tutorial notebook [`scDiscovery_cross_cell_state.ipynb`](scDiscovery_cross_cell_state.ipynb).

## Output

Running scDiscovery will generate the following files:

```text
output/
|   └── embedding.h5ad
figures/
|   ├── UMAP.png
|   ├── ROC.png
|   └── Density.png
```

The main output files include:

- `embedding.h5ad`: AnnData file containing learned cell embeddings
- `figures/`: visualization results

## Citation

If you find our code or method useful, please consider citing our work:

```bibtex
@article{huscDiscovery,
  title={scDiscovery: Entropy-guided adaptive delineation of novel cell types across single-cell omics},
  author={Zhaohui Hu, Chuan Li, Jialu Zhou, Siqi Tian, Yunhao Bai, Yu Zhang, Kunyao Zhu, Qun Yang, Nan Yin, Huibin Tan, Anxin Gu, Zhijiang Wan, Kaixiang Yang, Jinlong Shi, Guoli Yang, Xiaochun Cao, Liang Yang, Wenjing Yang, Long Lan, Yuansong Zeng, Li Shen, Mengzhu Wang},
  journal={},
  year={2026}
}
```
