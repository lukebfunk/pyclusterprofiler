import sharepathway
from scipy import stats
import numpy as np
import pandas as pd
import re
import warnings

from .utilities import map_database_gene_ids, correct_p

__all__ = [
    'compare_clusters'
    ]

def compare_clusters(df,grouping,enrichment_threshold=1,correction='fdr_bh',
    organism='hsa',database='KEGG',pathway_filters=None,verbose=True):

    if database=='KEGG':
        from .KEGG_utilities import (
            get_gene_id_mapping, 
            get_pathway_mapping, 
            get_gene_pathway_mapping
            )
    else:
        # todo: implement GO interface; probably can borrow functions from goatools (https://github.com/tanghaibao/goatools)
        raise NotImplementedError('Only `database`="KEGG" is currently implemented.')

    if pathway_filters is None:
        pathway_filter = lambda x: x
    elif callable(pathway_filters):
        pathway_filter = pathway_filters
    else:
        patterns = []
        if (database=='KEGG')&(organism=='hsa')&('remove_diseases' in pathway_filters):
            patterns.extend([r'.*hsa05.*',r'.*hsa0493[01234].*',r'.*hsa049[45]0.*',r'.*hsa015.*'])
        if (database=='KEGG')&(organism=='hsa')&('remove_organismal' in pathway_filters):
            patterns.extend([r'.*hsa046[12457].*',r'.*hsa04062.*',r'.*hsa049[1267].*',r'.*hsa04935.*',
                r'.*hsa03320.*',r'.*hsa042[67].*',r'.*hsa047.*',r'.*hsa043[268].*',r'.*hsa0421[123].*'])

        remaining = [f for f in pathway_filters if f not in ['remove_diseases','remove_organismal']]

        if len(remaining)>0:
            patterns.extend([p for p in remaining if isinstance(p,str)])

        pathway_filter = lambda x: not any([re.match(p,x) for p in patterns])

    gene_id_mapping = get_gene_id_mapping(organism)
    pathway_id_mapping = get_pathway_mapping(organism,filter=pathway_filter)
    gene_pathway_mapping = get_gene_pathway_mapping(organism, pathway_filter=pathway_filter)

    GENES = set()
    genelists = dict()
    clusters = dict()

    for i,(cluster,df_cluster_genes) in enumerate(df.groupby(grouping,sort=True)):
        genelists[cluster] = map_database_gene_ids(df_cluster_genes['gene_id'].tolist(),gene_id_mapping,verbose=verbose)
        clusters[i] = cluster
        GENES |= set(genelists[cluster])
    GENES = list(GENES)

    gene_cluster_matrix = genes2mat(genelists,GENES)

    pathways, pathway_counts, gene_pathway_matrix = sharepathway.linkpath2mat.linkpath2mat(GENES, gene_pathway_mapping)

    pathway_cluster_matrix = sharepathway.enrichment.enrichment(gene_cluster_matrix, gene_pathway_matrix)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore",category=RuntimeWarning)
        enrichment = ((pathway_cluster_matrix/gene_cluster_matrix.sum(axis=0))/
            (gene_pathway_matrix.sum(axis=0)/gene_pathway_matrix.shape[0]).T[:,None])

    results = []
    for p,c in zip(*(enrichment>enrichment_threshold).nonzero()):
        c_mask = np.zeros(enrichment.shape[1],dtype=bool)
        p_mask = np.zeros(enrichment.shape[0],dtype=bool)
        c_mask[c] = True
        p_mask[p] = True
        
        cg_mask = gene_cluster_matrix[:,c_mask].max(axis=1).astype(bool)
        
        notc_notp = (~gene_pathway_matrix[np.ix_(~cg_mask,p_mask)].astype(bool)).sum()
        notc_p = (gene_pathway_matrix[~cg_mask,p_mask]).sum()
        c_notp = (~gene_pathway_matrix[np.ix_(cg_mask,p_mask)].astype(bool)).sum()
        c_p = (gene_pathway_matrix[cg_mask,p_mask]).sum()
        
        odds_ratio, pvalue = stats.fisher_exact(np.array([[notc_p,c_p],[notc_notp,c_notp]]),alternative='less')

        cluster_pathway_ratio = pathway_cluster_matrix[p,c]/gene_cluster_matrix[:,c].sum()

        # prepare the table
        pathway_id = pathways[p]
        pathway_name = pathway_id_mapping[pathway_id].split(' - ')[0].strip()

        genes = [str(g) for g,include in zip(GENES,(gene_pathway_matrix[:,p].astype(bool)&gene_cluster_matrix[:,c].astype(bool))) 
                    if include]

        genesid = '+'.join([g.split(':')[1] for g in genes])
        map_url = "http://www.kegg.jp/pathway/"+pathway_id.split(':')[1]+'+'+genesid

        results.append({
            'cluster':clusters[c],
            'pathway':pathway_name,
            'pathway_id':pathway_id,
            'background_pathway_genes':pathway_counts[p],
            'background_pathway_ratio':gene_pathway_matrix[:,p].sum()/gene_cluster_matrix.shape[0],
            'cluster_pathway_genes':pathway_cluster_matrix[p,c],
            'cluster_pathway_ratio':pathway_cluster_matrix[p,c]/gene_cluster_matrix[:,c].sum(),
            'genes':genes,
            'observed/expected':enrichment[p,c],
            'pvalue':pvalue,
            'odds_ratio':odds_ratio,
            'map_url':map_url,
                       })

    df_results = pd.DataFrame(results)

    df_results['corrected_pvalue'] = (df_results
                                        .groupby(grouping,sort=False,as_index=False,group_keys=False)
                                        .apply(correct_p,method=correction)
                                     )

    df_results['-log10_corrected_pvalue'] = -np.log10(df_results['corrected_pvalue'])
    df_results['log2_observed/expected'] = np.log2(df_results['observed/expected'])

    return df_results


def genes2mat(cluster_genelists,full_genelist):
    G = len(full_genelist)
    C = len(cluster_genelists)
    genes_matrix = np.zeros((G,C))
    for c,genelist in enumerate(cluster_genelists.values()):
        for gene in genelist:
            genes_matrix[full_genelist.index(gene),c] = 1
    return genes_matrix