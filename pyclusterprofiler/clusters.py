import sharepathway
from scipy import stats
import numpy as np
import pandas as pd
import warnings
from functools import partial

from .utilities import map_database_gene_ids, correct_p

__all__ = [
    'compare_clusters'
    ]

def compare_clusters(df,grouping,correction='fdr_bh',
    organism='hsa',database='KEGG',exclude=None,force=False,verbose=True):

    if database=='KEGG':
        from .KEGG_utilities import (
            get_gene_id_mapping, 
            get_pathway_mapping, 
            get_gene_pathway_mapping
            )

        if organism==9606:
            organism='hsa'

    elif 'GO' in database:
        from .GO_utilities import (
            get_gene_id_mapping,
            get_pathway_mapping,
            get_gene_pathway_mapping
            )

        if organism=='hsa':
            organism=9606

        # uses GO-slim if 'slim' in `ontology`
        get_pathway_mapping = partial(get_pathway_mapping,ontology=database)

    else:
        raise NotImplementedError('Only `database`="KEGG" or "GO" is currently implemented.')

    gene_id_mapping = get_gene_id_mapping(organism, force=force)
    pathway_id_mapping = get_pathway_mapping(organism, exclude=exclude, force=force)

    if 'GO' in database:
        # when using GO, gene_pathway_mapping uses same download as gene_id_mapping
        force=False
    gene_pathway_mapping = get_gene_pathway_mapping(organism, annotations=list(pathway_id_mapping.keys()), force=force)

    GENES = set()
    genelists = dict()
    clusters = dict()

    df = df.astype({'gene_id':'str'})

    for i,(cluster,df_cluster_genes) in enumerate(df.groupby(grouping,sort=True)):
        genelists[cluster] = map_database_gene_ids(df_cluster_genes['gene_id'].tolist(),gene_id_mapping,verbose=verbose)
        clusters[i] = cluster
        GENES |= set(genelists[cluster])
    GENES = list(GENES)

    print(f'Analyzing {len(GENES)} genes, {df.groupby(grouping).ngroups} groups, {len(pathway_id_mapping)} pathways')

    gene_cluster_matrix = genes2mat(genelists,GENES)

    pathways, pathway_counts, gene_pathway_matrix = sharepathway.linkpath2mat.linkpath2mat(GENES, gene_pathway_mapping)

    pathway_cluster_matrix = sharepathway.enrichment.enrichment(gene_cluster_matrix, gene_pathway_matrix)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore",category=RuntimeWarning)
        enrichment = ((pathway_cluster_matrix/gene_cluster_matrix.sum(axis=0))/
            (gene_pathway_matrix.sum(axis=0)/gene_pathway_matrix.shape[0]).T[:,None])

    results = []
    # for p,c in zip(*(enrichment>enrichment_threshold).nonzero()):
    for p,c in zip(*np.indices(enrichment.shape).reshape(2,-1)):
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
        pathway_name = pathway_id_mapping[pathway_id]

        if database=="KEGG":
            pathway_name = pathway_name.split(' - ')[0].strip()

        genes = [str(g) for g,include in zip(GENES,(gene_pathway_matrix[:,p].astype(bool)&gene_cluster_matrix[:,c].astype(bool))) 
                    if include]

        result = {
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
            'odds_ratio':odds_ratio}

        if database=="KEGG":
            genesid = '+'.join([g.split(':')[1] for g in genes])
            map_url = "http://www.kegg.jp/pathway/"+pathway_id.split(':')[1]+'+'+genesid
            result['map_url'] = map_url

        results.append(result)

    if len(results) == 0:
        raise ValueError(f'No enriched pathways identified')

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