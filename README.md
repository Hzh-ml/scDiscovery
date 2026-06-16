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
- [Contact](#contact)

## Datasets

We provide easy access to the datasets used in this project.

The datasets will be released together with the project resources. After downloading the datasets, please organize the files as follows:

```text
scDiscovery/
├── data/
│   ├── dataset1
│   └── dataset2
|   └── ...
├── figures
├── scDiscovery.ipynb
├── README.md
└── requirements.txt
```

The input data should be stored in `.h5ad` format. Each `.h5ad` file should contain cells as observations and genes or features as variables.

A typical AnnData object should include:

- `adata.X`: expression or feature matrix
- `adata.obs`: cell-level metadata
- `adata.var`: gene-level or feature-level metadata
- `adata.obs["Batch"]`: batch or dataset information, if available
- `adata.obs["CellType"]`: cell type label, if available

For best performance, we recommend filtering low-quality cells and lowly detected genes or features before running scDiscovery.

## Installation

To reproduce scDiscovery, we suggest first creating a conda environment:

```bash
conda create -n scDiscovery python=3.8
conda activate scDiscovery
```

Then clone this repository:

```bash
git clone https://github.com/cquzys/scDiscovery.git
cd scDiscovery
```

Install the required packages:

```bash
pip install -r requirements.txt
```

If you use GPU acceleration, please install the PyTorch version that matches your CUDA version.

For example:

```bash
pip install torch torchvision torchaudio
```

## Usage

### Data preprocessing

Before running scDiscovery, please prepare the input single-cell dataset in `.h5ad` format.

A typical preprocessing command is:

```bash
python preprocess.py \
    --input data/raw/example.h5ad \
    --output data/processed/example_processed.h5ad
```

The preprocessing step may include quality control, normalization, feature selection, scaling, and construction of the processed AnnData object.

### Model training

After preprocessing, run the following command to train scDiscovery:

```bash
python train.py \
    --data data/processed/example_processed.h5ad \
    --save_path output/example
```

The trained model and intermediate results will be saved in the specified output directory.

### Representation extraction

To extract cell embeddings from a trained model, run:

```bash
python extract_embedding.py \
    --data data/processed/example_processed.h5ad \
    --model output/example/model.pt \
    --output output/example/embedding.h5ad
```

The output `.h5ad` file stores the learned cell representations and can be used for visualization and downstream analysis.

### Prediction

To apply a trained scDiscovery model to query data, run:

```bash
python predict.py \
    --model output/example/model.pt \
    --query data/processed/query_processed.h5ad \
    --output output/example_prediction
```

The prediction results will be saved in the output directory.

## Tutorial

### Tutorial 1: Basic single-cell discovery workflow

1. Install the required environment according to [Installation](#installation).
2. Prepare the input dataset in `.h5ad` format.
3. Place the dataset under the `data/raw/` folder.
4. Run preprocessing.
5. Train scDiscovery.
6. Extract cell embeddings.
7. Perform visualization and downstream analysis.

Example workflow:

```bash
conda activate scDiscovery

python preprocess.py \
    --input data/raw/example.h5ad \
    --output data/processed/example_processed.h5ad

python train.py \
    --data data/processed/example_processed.h5ad \
    --save_path output/example

python extract_embedding.py \
    --data data/processed/example_processed.h5ad \
    --model output/example/model.pt \
    --output output/example/embedding.h5ad

python predict.py \
    --model output/example/model.pt \
    --query data/processed/query_processed.h5ad \
    --output output/example_prediction
```

### Tutorial 2: Visualization of learned embeddings

After obtaining the learned embeddings, users can visualize the results using UMAP or t-SNE.

```bash
python visualize.py \
    --input output/example/embedding.h5ad \
    --output output/example/figures
```

The generated figures will be saved under:

```text
output/example/figures/
```

### Tutorial 3: Downstream analysis

The learned representations from scDiscovery can be used for downstream tasks, including:

- cell clustering
- cell type annotation
- batch effect analysis
- marker gene discovery
- trajectory or state discovery
- query-to-reference mapping

Users can customize downstream analysis based on the metadata stored in the AnnData object.

## Output

Running scDiscovery will generate the following files:

```text
output/example/
├── model.pt
├── embedding.h5ad
├── prediction.csv
├── train.log
└── figures/
```

The main output files include:

- `model.pt`: trained model weights
- `embedding.h5ad`: AnnData file containing learned cell embeddings
- `prediction.csv`: prediction or annotation results
- `train.log`: training log
- `figures/`: visualization results

## Citation

If you find our code or method useful, please consider citing our work:

```bibtex
@article{scDiscovery,
  title={scDiscovery: A Computational Framework for Single-Cell Discovery},
  author={Zhaohui Hu},
  journal={},
  year={2026}
}
```

## Contact

For questions, suggestions, or bug reports, please open an issue in this repository.

Maintainer: Zhaohui Hu
