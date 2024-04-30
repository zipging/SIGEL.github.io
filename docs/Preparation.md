

# Preparation

## Installation

Downloading SpaCEX code from https://github.com/WLatSunLab/SpaCEX

```pytho
git clone https://github.com/WLatSunLab/SpaCEX.git
```

Rename SpaCEX-main as SpaCEX.

## Download Data

 You can get the mouse hippocampus dataset [ssq-mHippo](https://singlecell.broadinstitute.org/single_cell/study/SCP815/sensitive-spatial-genome-wide-expression-profiling-at-cellular-resolution#study-summary. ), the human [dorsolateral prefrontal cortex](https://www.sciencedirect.com/topics/psychology/dorsolateral-prefrontal-cortex) datasets [10x-hDLPFC]( http://spatial.libd.org/spatialLIBD) are available. The human breast cancer dataset [10x-hBC]( https://support.10xgenomics.com/spatial-gene-expression/datasets/1.1.0/V1_Breast_Cancer_Block_A_Section_1) can be obtained. The mouse embryo dataset based on 10x Visium [10x-mEmb]( https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178636) can be found. The mouse embryo dataset based on SeqFISH [sqf-mEmb](https://crukci.shinyapps.io/SpatialMouseAtlas/) is obtainable. 

You can also obtain the example data [here](https://drive.google.com/drive/folders/1C3Gk-HVYp2dQh4id8H68M9p8IWEOIut_?usp=drive_link) required for SpaCEX.and make sure these data are organized in the following structure:
```
-SpaCEX
  -data
    -151676_10xvisium.h5ad
    -DLPFC_matrix_151676.dat
    -mEmb
      -10x_mEmb_matrix.dat
      -sqf_mEmb_adata.h5ad
      -sqf_mEmb_matrix.dat
  -model_pretrained
  -...
```

## Data in spacific task

### SpaCEX-ETC

```python
from SpaCEX.src.main.SpaCEX import SpaCEX
## get data on 10x and sqf
adata = SpaCEX_ETC.get_data(data='sqf', data_type='adata')
adata, key_m, dataset_m = SpaCEX_ETC.data_process(adata)
# key_m, dataset_m = SpaCEX_ETC.get_data(data='sqf', data_type='image')
key_v, dataset_v = SpaCEX_ETC.get_data(data='10x', data_type='image')
```

### SpaCEX-SVG

```python
## get adata and image data
adata= SpaCEX.get_data(sample_id='151676', data_type='adata')
dataset, adata = SpaCEX.data_process(adata)
```

### SpaCEX-SC

```python
## get adata and image data
adata= SpaCEX.get_data(sample_id='151676', data_type='adata')
dataset, adata = SpaCEX.data_process(adata)
```





