# 0 Welcome to SIGEL
We develop the **S**patially **I**nformed **G**ene **E**mbedding **L**earning (**SIGEL**) that can simultaneously identify spatially co-expressed genes and learn semantically meaningful gene embeddings from SRT data through a pretext task of gene clustering. SIGEL first employs an image encoder to transform the spatial expression maps of genes into gene embeddings modeled by a Student’s t mixture distribution (SMM). Subsequently, a discriminatively boosted gene clustering algorithm is applied on the posterior soft assignments of genes to the mixture components, iteratively adapting the parameters of the encoder and the SMM.

<div style="text-align: center;">
    <img src="logo_name.png" style="max-width: 50%; height: auto;">
</div>


## 0.1 Introduction

The limited exploration into spatial gene co-expression within tissues has been a significant bottleneck in fully harnessing the spatial genomic context for more insightful gene representations. To bridge this gap, we introduce **SIGEL**, a novel few-shot, self-supervised learning model tailored for the genomic field.

![png](Figure 1.png)

As shown in the figure above , **SIGEL** generates semantically meaningful gene Representations (SGRs) by identifying spatial gene co-expression patterns. Cofunctional and enrichment analyses of SGRs endorse their utility as genomic contexts, validated through relational semantics and exploration of gene functional ontology. Three novel SGR-based methods are proposed for enhancing FISH-based spatial transcriptomics, detecting spatially variable genes, and spatial clustering. Extensive real data results affirm the superior performance of these methods, highlighting the utility of SGRs for downstream tasks.

# 1 Preparation

## 1.1 Installation

To use SIGEL, please download code from https://github.com/WLatSunLab/SIGEL or conduct code below:

```python
git clone https://github.com/WLatSunLab/SIGEL.git
```

## 1.2 How to get help for SIGEL

Questions about SIGEL's code or experiments can be asked in the GitHub issue section at https://github.com/WLatSunLab/SIGEL.

## 1.3 Example Datasets

You can access a variety of datasets for spatial gene expression analysis:

- **Mouse Hippocampus Dataset ([ssq-mHippo](https://singlecell.broadinstitute.org/single_cell/study/SCP815/sensitive-spatial-genome-wide-expression-profiling-at-cellular-resolution#study-summary))**: Available for detailed study.
- **Human Dorsolateral Prefrontal Cortex Datasets ([10x-hDLPFC](http://spatial.libd.org/spatialLIBD))**: Information and datasets are accessible here.
- **Human Breast Cancer Dataset ([10x-hBC](https://support.10xgenomics.com/spatial-gene-expression/datasets/1.1.0/V1_Breast_Cancer_Block_A_Section_1))**: Obtainable for research purposes.
- **Mouse Embryo Dataset based on 10x Visium ([10x-mEmb](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178636))**: Available for comprehensive genetic analysis.
- **Mouse Embryo Dataset based on SeqFISH ([sqf-mEmb](https://crukci.shinyapps.io/SpatialMouseAtlas/))**: Can be accessed for advanced spatial analysis.

Example data required for SIGEL is available [here](https://drive.google.com/drive/folders/1C3Gk-HVYp2dQh4id8H68M9p8IWEOIut_?usp=drive_link). Please ensure that these data are properly organized as followes:

```
 . <SIGEL>
        ├── ...
        ├── <data>
        │   ├── 151676_10xvisium.h5ad
        │   ├── DLPFC_matrix_151676.dat
        │   └── <mEmb>
        │       ├── 10x_mEmb_matrix.dat
        │       ├── sqf_mEmb_adata.h5ad
        │       └── qf_mEmb_matrix.dat
        ├── <model_pretrained>
        │   │
        └── ...

```


# 2 Train model to obtain SIGEL-Generated Embeddings (SGRs) via ST dataset

In this task, we use 10x-visium data adata to generate embeddings SGRs.

## 2.1 Data preprocessing

The original `adata` is processed and converted into an image dataset `dataset`, with the gene names also extracted as `gene_name`
```python
from SIGEL.src.main.SIGEL import SIGEL
import warnings
warnings.filterwarnings("ignore")
```

```python
## get adata and image data
adata= SIGEL.get_data(sample_id='151676', data_type='adata')
dataset, adata = SIGEL.data_process(adata)
gene_name = adata.var.index.values
```

```python
## get SGRs
z, model = SIGEL.train(dataset, pretrain=True)
```

**Everything else is prepared and ready for deployment. If you have any other needs, please directly contact Wenlin Li at the email: zipging@gmail.com. **