import os
import anndata as ad
import networkx as nx
import scanpy as sc
import scglue
from matplotlib import rcParams
from scipy.io import mmread
import csv
import numpy as np

from itertools import chain

import anndata as ad
import itertools
import networkx as nx
import pandas as pd
import scanpy as sc
import scglue
import seaborn as sns
from matplotlib import rcParams

def load_rna(tissue):
    metadata = pd.read_csv("data/raw_data/rna_atac/metadata.csv", index_col=0)
    filepath = metadata.loc[tissue, "data_path"]
    rna_counts = mmread(os.path.join(filepath, "rna_counts.mtx"))
    rna_counts = rna_counts.transpose().tocsr()
    rna = ad.AnnData(rna_counts)
    #cell names
    barcodes = pd.read_csv(os.path.join(filepath, "barcode.csv"))
    rna.obs_names = barcodes["x"]
    #gene names
    genes = pd.read_csv(os.path.join(filepath, "genes.csv"))
    rna.var_names = genes["x"]

    # get cell type annotation
    ct_annotation = pd.read_csv(os.path.join(filepath, "ct_annotation.csv"))
    rna.obs["celltype"] = pd.Categorical(ct_annotation["x"])
    #set counts layer
    rna.layers["counts"] = rna.X.copy()
    #set origin
    rna.obs["domain"] = "scRNA-seq"
    return rna

def load_rna_no_annotation(tissue):
    metadata = pd.read_csv("data/raw_data/rna_atac/metadata.csv", index_col=0)
    filepath = metadata.loc[tissue, "data_path"]
    rna_counts = mmread(os.path.join(filepath, "rna_counts.mtx"))
    rna_counts = rna_counts.transpose().tocsr()
    rna = ad.AnnData(rna_counts)
    #cell names
    barcodes = pd.read_csv(os.path.join(filepath, "barcode.csv"))
    rna.obs_names = barcodes["x"]
    #gene names
    genes = pd.read_csv(os.path.join(filepath, "genes.csv"))
    rna.var_names = genes["x"]

    # get cell type annotation
    # ct_annotation = pd.read_csv(os.path.join(filepath, "ct_annotation.csv"))
    # rna.obs["celltype"] = pd.Categorical(ct_annotation["x"])
    #set counts layer
    rna.layers["counts"] = rna.X.copy()
    #set origin
    rna.obs["domain"] = "scRNA-seq"
    return rna

def load_atac(tissue):
    metadata = pd.read_csv("data/raw_data/rna_atac/metadata.csv", index_col=0)
    filepath = metadata.loc[tissue, "data_path"]

    atac_counts = mmread(os.path.join(filepath, "atac_counts.mtx"))
    atac_counts = atac_counts.transpose().tocsr()
    atac = ad.AnnData(atac_counts)

    #cell names
    barcodes = pd.read_csv(os.path.join(filepath, "barcode.csv"))
    atac.obs_names = barcodes["x"]

    #granges names
    granges = pd.read_csv(os.path.join(filepath, "peaks.csv"))
    atac.var_names = granges["x"]

    #set origin
    atac.obs["domain"] = "scATAC-seq"
    return atac

# def load_activity(tissue):
#     metadata = pd.read_csv("data/raw_data/rna_atac/metadata.csv", index_col=0)
#     filepath = metadata.loc[tissue, "data_path"]
#     activity_counts = mmread(os.path.join(filepath, "activity_counts.mtx"))
#     activity_counts = activity_counts.transpose().tocsr()
#     activity = ad.AnnData(activity_counts)
#     #cell names
#     barcodes = pd.read_csv(os.path.join(filepath, "barcode.csv"))
#     atac.obs_names = barcodes["x"]
#     #gene names
#     genes = pd.read_csv(os.path.join(filepath, "genes.csv"))
#     rna.var_names = genes["x"]


def subsample(rna, atac, ct_annotation, tissue, sub_num):
    #subsample the population
    metadata = pd.read_csv("data/raw_data/rna_atac/metadata.csv", index_col=0)
    filepath = metadata.loc[tissue, "data_path"]
    
    # Read activity counts and identify non-empty columns
    activity_counts = mmread(os.path.join(filepath, "activity_counts.mtx"))
    # non_empty_cols = activity_counts.sum(axis=0).nonzero()[1]
    
    # Get valid cells
    # valid_cells = rna.obs_names[non_empty_cols]
    cells = rna.obs_names
    subsample_size = int(len(cells) * sub_num)
    
    # Perform subsampling
    subsample_cells = np.random.choice(cells, size=subsample_size, replace=False)
    
    # Create subsampled copies
    rna_sub = rna[subsample_cells].copy()
    atac_sub = atac[subsample_cells].copy()
    ct_annotation_sub = ct_annotation.loc[subsample_cells]
    
    return (rna_sub, atac_sub, ct_annotation_sub, subsample_size)

def subsample_no_annotation(rna, atac, tissue, sub_num):

    #subsample the population
    metadata = pd.read_csv("data/raw_data/rna_atac/metadata.csv", index_col=0)
    filepath = metadata.loc[tissue, "data_path"]
    
    # Read activity counts and identify non-empty columns
    activity_counts = mmread(os.path.join(filepath, "activity_counts.mtx"))
    # non_empty_cols = activity_counts.sum(axis=0).nonzero()[1]
    
    # Get valid cells
    # valid_cells = rna.obs_names[non_empty_cols]
    cells = rna.obs_names
    subsample_size = int(len(cells) * sub_num)
    
    # Perform subsampling
    subsample_cells = np.random.choice(cells, size=subsample_size, replace=False)
    
    # Create subsampled copies
    rna_sub = rna[subsample_cells].copy()
    atac_sub = atac[subsample_cells].copy()
    
    return (rna_sub, atac_sub, subsample_size)

def feature_selection(var_genes_ordered, tissue, cor_method, var_num, total_num):
    
    metadata = pd.read_csv("data/raw_data/rna_atac/metadata.csv", index_col=0)
    species = metadata.loc[tissue, "species"]
    corr = pd.read_csv(os.path.join("output/conserved_features", species, ".".join([cor_method, "csv"])), index_col=0)
    corr = corr.sort_values(by='pearson_correlation', ascending=False)
    selected_genes = list()
    cor_num = 0
    if var_num==total_num:
        return (var_genes_ordered, 'Inf')
    for gene in corr.index:
        cor_num += 1
        if gene in var_genes_ordered:
            selected_genes.append(gene)
        if len(selected_genes)==total_num:
            break
    # BORDER CASE: ran out of CLIC_genes before reaching total_num
    # fill from var_genes_ordered that never appeared in clic_genes, preserving var_genes_ordered order
    remainder = [g for g in var_genes_ordered if (g not in corr.index) and (g not in selected_genes)]

    num_to_fill = total_num - len(selected_genes)
    if num_to_fill > 0:
        selected_genes.extend(remainder[:num_to_fill])
    return (selected_genes, cor_num)

