# sc_average_overlap
A Python package for comparing clusters in single-cell RNA-seq data using the average overlap metric. This package is designed to work seamlessly with the [Scanpy][] toolkit.

[Scanpy]: https://github.com/scverse/scanpy

## Getting Started

We recommend a Python version of 3.9.0+ in order to use the package.

### Installation

Currently, the sc_average_overlap package may be downloaded from TestPyPI. 

```
pip install -i https://test.pypi.org/simple/ sc-average-overlap
```

## Instructions for Use

The functions in this package assume an AnnData object that contains group labels for the cells, and that differential expression analysis with `sc.tl.rank_genes_groups()` has been performed.

First you may import the sc_average_overlap package as follows:

```
import sc_average_overlap as ao
```

### *make_ao_dendrogram()*
```
def make_ao_dendrogram(
	adata: AnnData,
	groupby: str,
	linkage_method: str = 'complete',
	de_key: str = 'rank_genes_groups',
	key: str = None, 
	genes_to_filter: List = None
) 
```

Given an AnnData object *adata* as input, this will first calculate pairwise average overlap scores for the clusters defined by the *groupby* parameter. You may specify a list of curated genes with the *genes_to_filter* argument, in which case all average overlap scores are based on the rankings of the user-specified gene list for each cluster.

Once this function is called, dendrogram information is saved into the AnnData. You may use Scanpy's plotting functions and specify the use of the average overlap dendrogram saved in the AnnData object.

An example function call:
```
ao.make_ao_dendrogram(adata, groupby='leiden')
ao.make_ao_dendrogram(adata, groupby='leiden', genes_to_filter=marker_gene_set) # providing a list of genes in marker_gene_set
```


### *plot_ao_dendrogram()*
```
def plot_ao_dendrogram(
	adata: AnnData,
	key: str
)
```

Once `make_ao_dendrogram` has been called, you may plot the resulting cluster tree. *key* should be the name of the observations grouping used during DE gene analysis and when calling `make_ao_dendrogram` to make the tree.

An example function call:
```
ao.plot_ao_dendrogram(adata, key='leiden')
```

### *plot_ao_heatmap()*
```
def plot_ao_heatmap(
	adata: AnnData,
	key: str,
	plot_zscores: bool = False,
	annot_decimal_format: str = '.1g'
)
```
You can plot a heatmap of all pairwise average overlap scores computed by calling the `plot_ao_heatmap` function. *key* should be the name of the observations grouping used during DE gene analysis and when calling `make_ao_dendrogram` to make the tree.

You can also plot the distances converted into z-scores. Average overlap follows a normal distribution when calculated on conjoint ranked lists - that is, the two ranked lists contain the same elements, just ranked differently. This gives a statistical interpretation of the resulting average overlap score. Only do this when specifying a specific marker gene set to base cluster rankings from, performed through providing a list of genes to the *genes_to_filter* argument when calling `make_ao_dendrogram`.

An example function call:
```
ao.plot_ao_heatmap(adata, key='leiden')
ao.plot_ao_heatmap(adata, key='leiden', plot_zscores=True)
```


## Authors

* **Christopher Thai** 

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

## Acknowledgments
* [Khiabanian Lab](https://khiabanian-lab.org)
* [Herranz Lab](http://www.herranzlab.org/)