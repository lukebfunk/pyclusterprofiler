import pandas as pd
from statsmodels.stats.multitest import multipletests

__all__ = [
    'map_database_gene_ids',
    'correct_p'
    ]

def map_database_gene_ids(ncbi_gene_ids,gene_id_mapping,verbose=True):
    gene_list = []

    for gene_id in ncbi_gene_ids:
        if gene_id in gene_id_mapping.keys():
            gene_list.append(gene_id_mapping[gene_id])
        else:
            if verbose:
                print(f'{gene_id} not found in database')

    return gene_list

def correct_p(df,pval_col='pvalue',method='fdr_bh'):
    return pd.Series(multipletests(df[pval_col],method=method)[1],
#                         columns=[f'corrected_{pval_col}'],
                        index=df.index
                       )