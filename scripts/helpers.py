"Helper functions for preprocessing"
import anndata
import numpy as np
import pandas as pd
from typing import Literal
from scipy.sparse import issparse


def normalize_smartseq(adata: anndata.AnnData, gene_len_file: str, column: Literal['mean', 'median'] = 'mean') -> anndata.AnnData:
    """Normalize SMART-seq sequencing by gene length

    Parameters
    ----------
    adata : anndata.AnnData
        Dataset
    gene_len : str
        Output file from gtftools
    column: Literal['mean', 'median'], default mean
        Column name which is used for gene length
    """
    print("SMART-SEQ: Normalization")

    gene_len = pd.read_table(gene_len_file).set_index('gene')

    common_genes = adata.var_names.intersection(gene_len.index)
    print(f"SMART-SEQ: Common genes {common_genes.shape[0]}")

    lengths = gene_len.loc[common_genes, column].values
    normalized = adata[:, common_genes].copy()
    if issparse(normalized.X):
        normalized.X = normalized.X.A
    
    normalized.X = normalized.X / lengths * np.median(lengths)
    normalized.X = np.rint(normalized.X)

    return normalized
