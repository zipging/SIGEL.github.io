## Step1: Load data

```python
from SpaCEX.src.main.SpaCEX import SpaCEX
from sklearn.preprocessing import MinMaxScaler
from scipy.cluster import hierarchy
import seaborn as sns
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
```

```python
## get adata and image data
adata= SpaCEX.get_data(sample_id='151676', data_type='adata')
dataset, adata = SpaCEX.data_process(adata)
gene_name = adata.var.index.values
```

## Step2: Generate SGEs

 ```python
 ## train model with pretrained model
 y_pred, SGEs, model = SpaCEX.train(dataset=dataset, pretrain=False)
 SGEs = SGEs.detach().cpu().numpy()
 ```

    adata2image: 100%|██████████| 18639/18639 [05:25<00:00, 57.27gene/s]

## Step3: Hierarchical clustering

```python
## select the gene family
key_select = ['KRT1', 'KRT5', 'KRT7', 'KRT86', 'KRT81', 'KRT83', 'KRT6B', 'KRT6A', 'KRT8', 'KRT23', 'KRT33B', 'KRT31', 'KRT37',
               'HLA-A', 'HLA-E', 'HLA-C', 'HLA-F', 'HLA-B',
              'HLA-DRA', 'HLA-DRB5', 'HLA-DRB1', 'HLA-DQA1', 'HLA-DQB1', 'HLA-DQA2', 'HLA-DQB2', 'HLA-DPA1', 'HLA-DPB1'
              ]
combined_data = zip(gene_name, SGEs)
emb = dict(combined_data)
emb = {key:emb[key] for key, value in emb.items() if key in key_select}
my_dict = {}
for key in emb.keys():
    my_dict[key] = np.array(emb[key]).flatten()
sorted_dict = dict(sorted(my_dict.items(), key=lambda x: key_select.index(x[0])))
sorted_dict_value = np.array(list(sorted_dict.values()))
scaler = MinMaxScaler()
```

```python
## norm the data
normalized_array = scaler.fit_transform(sorted_dict_value.T).T
my_dict_norm = {}
for i in range(len(sorted_dict)):
    key = list(sorted_dict.keys())[i]
    my_dict_norm[key] = normalized_array[i]
```

```python
## plot hierarchical clustering results
df = pd.DataFrame.from_dict(my_dict_norm, orient='index')

sns.set(font_scale=2.5, font='sans-serif')
row_colors = ['#00C5CD'] * 13 + ['#1E90FF'] * 5 + ['#CD3333'] * 9

## hierarchical clustering
row_linkage = hierarchy.linkage(df, method='average')
col_linkage = hierarchy.linkage(df.transpose(), method='average')

# plot heatmap
sns.set(font_scale=1., font='sans-serif' )
sns.clustermap(df, row_linkage=row_linkage, col_linkage=col_linkage, cmap='coolwarm', row_colors=[row_colors], figsize=(10,8))

plt.show()
```


![png](output_1_0.png)



