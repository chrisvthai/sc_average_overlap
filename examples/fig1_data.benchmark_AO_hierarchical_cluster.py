import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import argparse

from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score, fowlkes_mallows_score
from sklearn.metrics.cluster import pair_confusion_matrix
from sklearn.preprocessing import LabelEncoder
from scipy.cluster.hierarchy import fcluster
from scipy.spatial import distance

sys.path.insert(0, '../')
from tools import average_overlap as ao

# Get arguments first
parser = argparse.ArgumentParser()
parser.add_argument("--resolution", type=float)
parser.add_argument("--n_markers", type=int)
parser.add_argument("--n_iters", type=int)
args = parser.parse_args()

resolution = args.resolution
n_markers = args.n_markers
n_iters = args.n_iters


### Read and preprocess data
adata = sc.read_10x_mtx('ZhengMix8eq', prefix='ZhengMix8eq_')
df_celllabels = pd.read_csv('ZhengMix8eq/ZhengMix8eq_KnownCellLabels.tsv', sep='\t', header=None)
adata.obs['label'] = pd.Categorical(df_celllabels[0])

t_cell_labels = ['memory.t', 'naive.cytotoxic', 'naive.t', 'regulatory.t', 'cd4.t.helper']
adata = adata[adata.obs['label'].isin(t_cell_labels)]

# Standard preprocessing
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

#sc.pp.highly_variable_genes(adata, flavor='seurat_v3', n_top_genes=5000)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

adata.raw = adata

# PCA + UMAP
#sc.tl.pca(adata)
adata.obsm['X_pca'] = pd.read_csv('ZhengMix8eq_Piccolo_PCA.csv', index_col=None).to_numpy().astype(np.float32)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=50)

#### Main for loop
amis = []
aris = []
fms_arr = []
tps = []
tns = []
fps = []
fns = []

labels_dict = {'orig_labels': adata.obs['label'].tolist()}
distances_dict = {}

for i in range(0, n_iters):
    # 1. Cluster with Leiden
    sc.tl.leiden(adata, resolution=resolution, random_state=i)
    
    # 2. Generate AO dendrogram
    marker_gene_set = set()
    cat = 'leiden'
    sc.tl.rank_genes_groups(adata, groupby=cat, method='wilcoxon', tie_correct=True)

    for c in adata.obs[cat].cat.categories:
        marker_gene_set.update(list(ao.get_cluster_markers(adata, cluster_label=c, n_genes=n_markers)))

    ao.make_ao_dendrogram(adata, groupby=cat, genes_to_filter=marker_gene_set)
    

    # 3. Automatically merge hierarchical clusters into 5 final clusters
    linkage = adata.uns['dendrogram_ao_leiden']['linkage']
    merged_cluster_labels = fcluster(linkage, 5, criterion='maxclust').astype('str') # Merge to 5 clusters

    merge_dict = {}
    for x, y in zip(adata.uns['dendrogram_ao_leiden']['categories'], merged_cluster_labels):
        merge_dict[x] = y

    adata.obs['leiden_merge'] = adata.obs['leiden']
    adata.obs['leiden_merge'] = (
        adata.obs['leiden_merge']
        .map(lambda x: merge_dict.get(x, x)) # get() with 2 arguments returns the key itself if it doesn't exist in dict
        .astype('category')
    )

    # 4. Calculate AMI and ARI 
    # AMI + ARI calculation
    cols = ['label', 'leiden_merge', ]

    labels = adata.obs[cols]
    le = LabelEncoder()
    for c in cols:
        labels['{}_encoded'.format(c)] = le.fit_transform(labels[c])

    ari = adjusted_rand_score(labels['label_encoded'], labels['leiden_merge_encoded'])
    ami = adjusted_mutual_info_score(labels['label_encoded'], labels['leiden_merge_encoded'])
    fms = fowlkes_mallows_score(labels['label_encoded'], labels['leiden_merge_encoded'])
    
    conf_mat = pair_confusion_matrix(labels['label_encoded'], labels['leiden_merge_encoded'])
    tp = conf_mat[0, 0]
    tn = conf_mat[1, 1]
    fp = conf_mat[0, 1]
    fn = conf_mat[1, 0]
    
    aris.append(ari)
    amis.append(ami)
    fms_arr.append(fms)
    tps.append(tp)
    tns.append(tn)
    fps.append(fp)
    fns.append(fn)

    # 5.  Recompute the tree with the same parameters, except with the merged clusters, and get pairwise distances
    marker_gene_set = set()
    cat = 'leiden_merge'
    sc.tl.rank_genes_groups(adata, groupby=cat, method='wilcoxon', tie_correct=True)

    for c in adata.obs[cat].cat.categories:
        marker_gene_set.update(list(ao.get_cluster_markers(adata, cluster_label=c, n_genes=n_markers)))

    ao.make_ao_dendrogram(adata, groupby=cat, genes_to_filter=marker_gene_set)

    distances = adata.uns['ao_leiden_merge']['ao_heatmap']
    distances = distance.squareform(1 - distances)
    distances_dict[i] = distances.tolist()
    
    # 6.  Save a table with the 1st column being ground truth, and the nth column being the merged labels
    labels_dict[str(i)] = adata.obs['leiden_merge'].tolist()

pd.DataFrame({'ARI': aris, 'AMI': amis, 'FMS': fms_arr, 'TP': tps, 'TN': tns, 'FP': fps, 'FN': fns}).to_csv('ao_dendro_benchmarks_r={}_nmarkers={}_SCORES.csv'.format(resolution, n_markers), index=False)
pd.DataFrame.from_dict(distances_dict, orient='index').to_csv('ao_dendro_benchmarks_r={}_nmarkers={}_DISTANCES.csv'.format(resolution, n_markers), index=False)
pd.DataFrame.from_dict(labels_dict).to_csv('ao_dendro_benchmarks_r={}_nmarkers={}_LABELS.csv'.format(resolution, n_markers), index=False)

