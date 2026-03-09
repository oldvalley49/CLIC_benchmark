import anndata as ad
import networkx as nx
import scanpy as sc
import scglue
from matplotlib import rcParams
from scipy.io import mmread
from scipy.io import mmwrite
from scipy import sparse
import csv
import numpy as np
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neighbors import KNeighborsRegressor
from itertools import chain

import pandas as pd
import scglue
import seaborn as sns
import sys
import os
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix
from utility import *
scglue.plot.set_publication_params()

# command-line arguments
args = sys.argv[1:]
tissue = args[0]
var_num = int(args[1])
total_num = int(args[2])
sub_num = float(args[3])
cor_method = args[4]
index = int(args[5])
activity_model = args[6]

data_dir = "data/processed_data/rna_atac/"+tissue
rna = load_rna(tissue)
atac = load_atac(tissue)
metadata = pd.read_csv("data/raw_data/rna_atac/metadata.csv", index_col=0)
annotation = metadata.loc[tissue, "annotation"]
species = metadata.loc[tissue, "species"]
gtf_fp = os.path.join("data/raw_data/gtf", ".".join([annotation, "gtf"]))

#scglue.config.BEDTOOLS_PATH = "/opt/homebrew/bin"

# load cell type annotation
ct_annotation = pd.read_csv(os.path.join(data_dir, "ct_annotation.csv"))
barcodes = pd.read_csv(os.path.join(data_dir, "barcode.csv"), index_col=0)
barcodes = barcodes['x']
barcodes.reset_index(drop=True, inplace=True)
ct_annotation['barcodes'] = barcodes
ct_annotation.set_index('barcodes', inplace=True)
#subsample the population

# cells = rna.obs_names
# subsample_size = int(len(cells)*sub_num)
# subsample_cells = np.random.choice(cells, size=subsample_size, replace=False)
# rna = rna[subsample_cells].copy()
# atac = atac[subsample_cells].copy()
# ct_annotation = ct_annotation.loc[subsample_cells]
rna,atac,ct_annotation,subsample_size = subsample(rna, atac, ct_annotation, tissue, sub_num)
# exit if there are more more var_num than the number of features
if var_num > rna.n_vars:
    raise ValueError("there are more requested intial_var_features than there are features in the original dataset. ")
# RNA preprocessing
sc.pp.highly_variable_genes(rna, n_top_genes=var_num, flavor="seurat_v3")

var_genes = rna.var['variances_norm'].sort_values(ascending=False).index[0:var_num]
subset_genes,cor_num = feature_selection(var_genes, tissue, cor_method, var_num, total_num)

# generate run_id
rna.var['highly_variable'] = rna.var.index.isin(subset_genes)
test = rna.var[rna.var['highly_variable'] == True]
print(test)
run_id = "_".join(["glue", cor_method, str(var_num), str(cor_num), str(subsample_size), activity_model, str(index)])
print(run_id)

sc.pp.normalize_total(rna)
sc.pp.log1p(rna)
norm_data = rna.X.copy()
sc.pp.scale(rna)
sc.tl.pca(rna, n_comps=100, svd_solver="auto")
sc.pp.neighbors(rna, metric="cosine")
scglue.data.get_gene_annotation(
    rna, gtf=gtf_fp,
    gtf_by="gene_name"
)
valid_chroms = rna.var['chrom'].dropna().index
rna = rna[:, valid_chroms]
# ATAC preprocessing
scglue.data.lsi(atac, n_components=100, n_iter=15)
sc.pp.neighbors(atac, use_rep="X_lsi", metric="cosine")
split = atac.var_names.str.split(r"[-]")
atac.var["chrom"] = split.map(lambda x: x[0])
atac.var["chromStart"] = split.map(lambda x: x[1]).astype(int)
atac.var["chromEnd"] = split.map(lambda x: x[2]).astype(int)


guidance = scglue.genomics.rna_anchored_guidance_graph(rna, atac)
scglue.graph.check_graph(guidance, [rna, atac])
scglue.models.configure_dataset(
    rna, "NB", use_highly_variable=True,
    use_layer="counts", use_rep="X_pca"
)
scglue.models.configure_dataset(
    atac, "NB", use_highly_variable=True,
    use_rep="X_lsi"
)
guidance_hvf = guidance.subgraph(chain(
    rna.var.query("highly_variable").index,
    atac.var.query("highly_variable").index
)).copy()
progress_dir = os.path.join("glue", run_id)
glue = scglue.models.fit_SCGLUE(
    {"rna": rna, "atac": atac}, guidance_hvf,
    fit_kws={"directory": progress_dir}
)
# model_filename = run_id + ".dill"
# glue.save(os.path.join("models/glue", model_filename))
dx = scglue.models.integration_consistency(
    glue, {"rna": rna, "atac": atac}, guidance_hvf
)
rna.obsm["X_glue"] = glue.encode_data("rna", rna)
atac.obsm["X_glue"] = glue.encode_data("atac", atac)
combined = ad.concat([rna, atac])
sc.pp.neighbors(combined, use_rep="X_glue", metric="cosine")
sc.tl.umap(combined)
fig = sc.pl.umap(combined, color=["domain"], return_fig=True)
directory_path = os.path.join("plots/paired/rna_atac", tissue, "umap")
os.makedirs(directory_path, exist_ok=True)
plt.savefig(os.path.join(directory_path, ".".join([run_id, "jpeg"])))
feature_embeddings = glue.encode_graph(guidance_hvf)
feature_embeddings = pd.DataFrame(feature_embeddings, index=glue.vertices)
feature_embeddings.iloc[:5, :]
rna.varm["X_glue"] = feature_embeddings.reindex(rna.var_names).to_numpy()
atac.varm["X_glue"] = feature_embeddings.reindex(atac.var_names).to_numpy()
coembedding = pd.DataFrame(combined.obsm["X_glue"])
coembedding.columns = ["PC" + str(i) for i in range(1, 51)]
coembedding.index = combined.obs_names
directory_path = os.path.join("output/paired/rna_atac", tissue, "coembed")
os.makedirs(directory_path, exist_ok=True)
# save coembedding
coembedding.to_csv(os.path.join(directory_path, ".".join([run_id, "csv"])))

# transfer annotations using coembedding space with KNN classifier

# load annotation information
rna_embed = coembedding.iloc[:subsample_size, :]
atac_embed = coembedding.iloc[subsample_size:, :]
ct_truth = ct_annotation['x']
transfer_annotation = KNeighborsClassifier(n_neighbors=30)
transfer_annotation.fit(rna_embed, ct_truth)
ct_predicted = transfer_annotation.predict(atac_embed)

celltypes = pd.DataFrame({
    'barcodes': ct_truth.index,
    'ct_truth': ct_truth.values,
    'predicted.id': ct_predicted
})
directory_path = os.path.join("output/paired/rna_atac", tissue, "ct_annotation")
os.makedirs(directory_path, exist_ok=True)
celltypes.to_csv(os.path.join(directory_path, ".".join([run_id, "csv"])))

conf_matrix = confusion_matrix(celltypes['ct_truth'], celltypes['predicted.id'], labels=np.unique(celltypes['ct_truth']))
conf_matrix_normalized = conf_matrix / conf_matrix.sum(axis=1, keepdims=True)

conf_matrix_df = pd.DataFrame(
    conf_matrix_normalized, 
    index=np.unique(celltypes['ct_truth']), 
    columns=np.unique(celltypes['ct_truth'])
)

# conf_matrix_df = pd.DataFrame(
#     conf_matrix, 
#     index=np.unique(celltypes['ct_truth']), 
#     columns=np.unique(celltypes['ct_truth'])
# )
plt.figure(figsize=(12, 10))
sns.heatmap(conf_matrix_df, annot=False, cmap='Blues', fmt='g')
directory_path = os.path.join("plots/paired/rna_atac", tissue, "annotation_accuracy")
os.makedirs(directory_path, exist_ok=True)
plt.savefig(os.path.join(directory_path,".".join([run_id, "jpeg"])), dpi=300)
