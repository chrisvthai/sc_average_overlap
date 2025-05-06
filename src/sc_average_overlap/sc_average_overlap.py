import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import random

from scipy.spatial import distance
from scipy.cluster import hierarchy
from anndata import AnnData
from typing import List, Optional, Union


def _compute_overlap(l1: Union[List, np.ndarray], l2: Union[List, np.ndarray]):
	assert len(l1) == len(l2), "lists must be equal in length"
	return len(set(l1).intersection(l2)) / len(l1)

def _average_overlap(l1: Union[List, np.ndarray], l2: Union[List, np.ndarray]):
	assert type(l1) in [list, np.ndarray], "first array must be np.ndarray or List"
	assert type(l2) in [list, np.ndarray], "second array must be np.ndarray or List"
	assert len(l1) == len(l2), "lists must be equal in length"

	total_depth = len(l1)
	ao = 0

	for d in range(1, total_depth+1):
		o = _compute_overlap(l1[:d], l2[:d])
		ao += o

	ao /= total_depth
	return ao


def get_cluster_markers(
	adata: AnnData,
	cluster_label: str,
	genes_to_filter: List = None,
	n_genes: int = 50,
	key: str = 'rank_genes_groups'
):
	"""
	A helper function for extracting marker genes for a given cluster

	Given sc.tl.rank_genes_groups() has already been run for one-vs-rest differential expression,
	retrieve markers of the given cluster.

	If a list of genes is specified, then this gives a list of those genes ranked by their differential expression for that cluster
	"""
	all_markers = adata.uns[key]['names'][cluster_label]

	if genes_to_filter is not None:
		all_markers = [g for g in all_markers if g in genes_to_filter]
		n_genes = len(all_markers)
	return all_markers[:n_genes]

def make_ao_dendrogram(
	adata: AnnData,
	groupby: str,
	linkage_method: str = 'complete',
	de_key: str = 'rank_genes_groups',
	key: str = None, # user can specify a specific key to save AO data into. Otherwise groupby name is used
	genes_to_filter: List = None,
):

	cluster_labels = adata.obs[groupby].cat.categories

	# Distance and z-score heatmaps
	ao_heatmap = np.empty((len(cluster_labels), len(cluster_labels)), dtype=np.double)
	for x, k in enumerate(cluster_labels):
		for y, j in enumerate(cluster_labels):
			l1 = get_cluster_markers(adata, k, key=de_key, genes_to_filter=genes_to_filter)
			l2 = get_cluster_markers(adata, j, key=de_key, genes_to_filter=genes_to_filter)
			
			ao_heatmap[x, y] = _average_overlap(l1, l2)
	
	# Depending on the gene list given and its length, fit a normal distribution for AO on random shuffles
	def get_normal_parameters(
		gene_list: List,
		n_iter_random_samples: int = 500,
	):
		from scipy.stats import norm
		## sample AO from randomly shuffled lists, and fit a normal
		l1 = gene_list.copy()
		l2 = gene_list.copy()

		ao_score_dist = []
		for i in range(0, n_iter_random_samples):
			random.shuffle(l1)
			random.shuffle(l2)
			ao = _average_overlap(l1, l2)
			ao_score_dist.append(ao)

		mu, std = norm.fit(ao_score_dist)
		return mu, std
	
	mu, std = get_normal_parameters(l1)
	make_zscore = lambda x: (x - mu) / std
	zscore_heatmap = make_zscore(ao_heatmap)

	dat_ao = dict(
		ao_heatmap = ao_heatmap,
		zscore_heatmap = zscore_heatmap,
		mu = mu,
		std = std,
		cluster_labels = cluster_labels
	)

	if key is None:
		ao_key = 'ao_{}'.format(groupby)
	else:
		ao_key = 'ao_{}'.format(key)
	adata.uns[ao_key] = dat_ao

	# Hierarchical clustering
	dist_heatmap = 1.0 - ao_heatmap
	distances = distance.squareform(dist_heatmap)
	Z = hierarchy.linkage(distances, method=linkage_method, )
	dn = hierarchy.dendrogram(Z, labels=list(cluster_labels), no_plot=True)

	dat_dendro = dict(
		linkage = Z,
		groupby = [groupby],
		linkage_method = linkage_method,
		categories = list(cluster_labels),
		categories_ordered = dn["ivl"],
		categories_idx_ordered = dn["leaves"],
		dendrogram_info = dn,
		distance_matrix = dist_heatmap
	)

	if key is None:
		dendrogram_key = 'dendrogram_ao_{}'.format(groupby)
	else:
		dendrogram_key = 'dendrogram_ao_{}'.format(key)
	adata.uns[dendrogram_key] = dat_dendro


# Separate plotting function for dendrograms, outside of scanpy's
def plot_ao_dendrogram(
	adata: AnnData,
	key: str,
):
	ao_key = 'ao_{}'.format(key)
	dendrogram_key = 'dendrogram_ao_{}'.format(key)

	if dendrogram_key not in adata.uns or ao_key not in adata.uns: #TODO: replace print statement with an actual warning
		print(
			"Dendrogram or ao data not found (using key={}).\nRun average_overlap.make_ao_dendrogram() first."
			.format(key)
		)
		return

	dat_dendro = adata.uns[dendrogram_key]
	dat_ao = adata.uns[ao_key]
	mu = dat_ao['mu']
	std = dat_ao['std']

	fig, ax = plt.subplots()
	dn = hierarchy.dendrogram(
			dat_dendro['linkage'], 
			labels=dat_dendro['categories'], 
			ax=ax,
			color_threshold=0,
			above_threshold_color='black'
		)
	
	def ao_to_z(ao):
		return (1 - ao - mu) / std

	def z_to_ao(z):
		return 1 - (z * std) - mu

	ax.set_xlabel('Cluster label')
	ax.set_ylabel('AO Distance')
	secax = ax.secondary_yaxis('right', functions=(ao_to_z, z_to_ao))
	secax.set_ylabel('z-score')

	return fig, ax

def plot_ao_heatmap(
	adata: AnnData,
	key: str,
	plot_zscores: bool = False,
	annot_decimal_format: str = '.1g'
):
	ao_key = 'ao_{}'.format(key)
	dendrogram_key = 'dendrogram_ao_{}'.format(key)

	if ao_key not in adata.uns: #TODO: replace print statement with an actual warning
		print(
			"Ao data not found (using key={}).\nRun average_overlap.make_ao_dendrogram() first."
			.format(key)
		)
		return

	dat_ao = adata.uns[ao_key]
	dat_dendro = adata.uns[dendrogram_key]
	if plot_zscores:
		heatmap = dat_ao['zscore_heatmap']
	else:
		heatmap = dat_ao['ao_heatmap']

	linkage = dat_dendro['linkage']

	fig, ax = plt.subplots()

	ax = sns.heatmap(heatmap, annot=True, fmt=annot_decimal_format, annot_kws={"size": 10},
				xticklabels=dat_ao['cluster_labels'], yticklabels=dat_ao['cluster_labels']
	)

	return fig, ax
	