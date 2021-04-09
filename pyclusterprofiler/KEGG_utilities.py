import sharepathway

__all__ = [
	'get_gene_id_mapping',
	'get_pathway_mapping',
	'get_gene_pathway_mapping'
	]

def get_gene_id_mapping(organism='hsa'):
	return sharepathway.geneIDconv.geneIDconv(species=organism)

def get_pathway_mapping(organism='hsa'):
	pathway_tuples=sharepathway.parse_kegg.Request('list','pathway',organism)
	return {pathway_id:name for pathway_id,name in pathway_tuples}

def get_gene_pathway_mapping(organism='hsa'):
    return sharepathway.parse_kegg.Request('link','path',organism)