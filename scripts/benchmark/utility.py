import numpy as np
import pandas as pd
import os
import sys
import re
import scipy
import sklearn.metrics

from typing import Tuple

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.spatial
import sklearn.metrics
import sklearn.neighbors
from anndata import AnnData
from scipy.sparse.csgraph import connected_components
#from .typehint import RandomState
#from .utils import get_rs


def load_coembed(file_path):
    # load embeddings
    coembed = pd.read_csv(file_path, index_col=0)
    cell_num = int(coembed.shape[0]/2)
    # split into RNA and ATAC dataset
    rna = coembed[0:cell_num]
    atac = coembed[cell_num:]
    return rna, atac

def load_coembed_unpaired(file_path, rna_num):
    # load embeddings
    coembed = pd.read_csv(file_path, index_col=0)
    print(coembed)
    # split into RNA and ATAC dataset
    rna = coembed[0:rna_num]
    atac = coembed[rna_num:]
    return rna, atac
    
def load_rna(file_path):
    rna = pd.read_csv(file_path, index_col=0)
    cell_num = int(rna.shape[0]/2)
    return rna
# measures alignment quality
def foscttm(x, y):
    if x.shape != y.shape:
        raise ValueError("Shapes do not match!")
    d = scipy.spatial.distance_matrix(x, y)
    foscttm_x = (d < np.expand_dims(np.diag(d), axis=1)).mean(axis=1)
    foscttm_y = (d < np.expand_dims(np.diag(d), axis=0)).mean(axis=0)
    fracs = []
    for i in range(len(foscttm_x)):
        fracs.append((foscttm_x[i] + foscttm_y[i]) / 2)
    return 1-np.mean(fracs).round(4)
# measures cell type annotation accuracy
def ari(x, y):
    """
    x: true labels
    y: predicted labels
    """
    return sklearn.metrics.adjusted_rand_score(x, y)

# raw accuracy percentage
def accuracy(x, y):
    """
    x: true labels
    y: predicted labels
    """
    return sklearn.metrics.accuracy_score(x, y)
def mean_average_precision(
        x: np.ndarray, y: np.ndarray, neighbor_frac: float = 0.01, **kwargs
) -> float:
    r"""
    Mean average precision

    Parameters
    ----------
    x
        Coordinates
    y
        Cell type labels
    neighbor_frac
        Nearest neighbor fraction
    **kwargs
        Additional keyword arguments are passed to
        :class:`sklearn.neighbors.NearestNeighbors`

    Returns
    -------
    map
        Mean average precision
    """
    k = max(round(y.shape[0] * neighbor_frac), 1)
    nn = sklearn.neighbors.NearestNeighbors(
        n_neighbors=min(y.shape[0], k + 1), **kwargs
    ).fit(x)
    nni = nn.kneighbors(x, return_distance=False)
    match = np.equal(y[nni[:, 1:]], np.expand_dims(y, 1))
    return np.apply_along_axis(_average_precision, 1, match).mean().item()



def _average_precision(match: np.ndarray) -> float:
    if np.any(match):
        cummean = np.cumsum(match) / (np.arange(match.size) + 1)
        return cummean[match].mean().item()
    return 0.0


def normalized_mutual_info(x: np.ndarray, y: np.ndarray, **kwargs) -> float:
    r"""
    Normalized mutual information with true clustering

    Parameters
    ----------
    x
        Coordinates
    y
        Cell type labels
    **kwargs
        Additional keyword arguments are passed to
        :func:`sklearn.metrics.normalized_mutual_info_score`

    Returns
    -------
    nmi
        Normalized mutual information

    Note
    ----
    Follows the definition in `OpenProblems NeurIPS 2021 competition
    <https://openproblems.bio/neurips_docs/about_tasks/task3_joint_embedding/>`__
    """
    x = AnnData(X=x, dtype=x.dtype)
    sc.pp.neighbors(x, n_pcs=0, use_rep="X")
    nmi_list = []
    for res in (np.arange(20) + 1) / 10:
        sc.tl.leiden(x, resolution=res)
        leiden = x.obs["leiden"]
        nmi_list.append(sklearn.metrics.normalized_mutual_info_score(
            y, leiden, **kwargs
        ).item())
    return max(nmi_list)



def asw(x: np.ndarray, y: np.ndarray, **kwargs) -> float:
    r"""
    Cell type average silhouette width

    Parameters
    ----------
    x
        Coordinates
    y
        Cell type labels
    **kwargs
        Additional keyword arguments are passed to
        :func:`sklearn.metrics.silhouette_score`

    Returns
    -------
    asw
        Cell type average silhouette width

    Note
    ----
    Follows the definition in `OpenProblems NeurIPS 2021 competition
    <https://openproblems.bio/neurips_docs/about_tasks/task3_joint_embedding/>`__
    """
    return (sklearn.metrics.silhouette_score(x, y, **kwargs).item() + 1) / 2



def graph_connectivity(
        x: np.ndarray, y: np.ndarray, **kwargs
) -> float:
    r"""
    Graph connectivity

    Parameters
    ----------
    x
        Coordinates
    y
        Cell type labels
    **kwargs
        Additional keyword arguments are passed to
        :func:`scanpy.pp.neighbors`

    Returns
    -------
    conn
        Graph connectivity
    """
    x = AnnData(X=x, dtype=x.dtype)
    sc.pp.neighbors(x, n_pcs=0, use_rep="X", **kwargs)
    conns = []
    for y_ in np.unique(y):
        x_ = x[y == y_]
        _, c = connected_components(
            x_.obsp['connectivities'],
            connection='strong'
        )
        counts = pd.value_counts(c)
        conns.append(counts.max() / counts.sum())
    return np.mean(conns).item()

def asw_batch(
        x: np.ndarray, y: np.ndarray, ct: np.ndarray, **kwargs
) -> float:
    r"""
    Batch average silhouette width

    Parameters
    ----------
    x
        Coordinates
    y
        Batch labels
    ct
        Cell type labels
    **kwargs
        Additional keyword arguments are passed to
        :func:`sklearn.metrics.silhouette_samples`

    Returns
    -------
    asw_batch
        Batch average silhouette width

    Note
    ----
    Follows the definition in `OpenProblems NeurIPS 2021 competition
    <https://openproblems.bio/neurips_docs/about_tasks/task3_joint_embedding/>`__
    """
    s_per_ct = []
    for t in np.unique(ct):
        mask = ct == t
        try:
            s = sklearn.metrics.silhouette_samples(x[mask], y[mask], **kwargs)
        except ValueError:  # Too few samples
            s = 0
        s = (1 - np.fabs(s)).mean()
        s_per_ct.append(s)
    return np.mean(s_per_ct).item()

def knn_auc(F: np.ndarray, G: np.ndarray, k_percent: float = 0.01) -> float:
    """
    Computes the kNN-AUC score, measuring the average percentage overlap between
    the neighborhoods of scRNA-seq measurements (F) and scATAC-seq measurements (G).
    
    Parameters
    ----------
    F : np.ndarray
        scRNA-seq measurements (cells × features)
    G : np.ndarray
        scATAC-seq measurements (cells × features)
    k : int, optional
        Number of neighbors to consider as percentage fo the entire population, default is 1%.
    
    Returns
    -------
    float
        The kNN-AUC score (higher is better).
    """
    assert F.shape[0] == G.shape[0], "F and G must have the same number of cells"
    n_cells = F.shape[0]
    n_neighbors  = int(n_cells*k_percent)
    # Compute kNN in F
    nn_F = sklearn.neighbors.NearestNeighbors(n_neighbors=n_neighbors+1).fit(F)
    neighbors_F = nn_F.kneighbors(F, return_distance=False)[:, 1:]  # Exclude self
    
    # Compute kNN in G
    nn_G = sklearn.neighbors.NearestNeighbors(n_neighbors=n_neighbors+1).fit(G)
    neighbors_G = nn_G.kneighbors(G, return_distance=False)[:, 1:]  # Exclude self
    
    # Compute overlap for each cell
    overlaps = [len(set(neighbors_F[i]) & set(neighbors_G[i])) / n_neighbors for i in range(n_cells)]
    
    # Return average percentage overlap
    return np.mean(overlaps)

