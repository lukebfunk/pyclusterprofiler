import os
import goatools.base
from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.obo_parser import GODag

__all__ = [
	'get_gene_id_mapping',
	'get_pathway_mapping',
	'get_gene_pathway_mapping'
	]

def get_gene_id_mapping(organism=9606, force=False):
	if force & (os.path.isfile('gene2go')):
		os.remove('gene2go')

	gene2go = Gene2GoReader(goatools.base.download_ncbi_associations(), taxids=[organism]).get_id2gos_nss()

	return {str(gene):gene for gene in gene2go.keys()}

def get_pathway_mapping(organism=9606, ontology='basic', exclude=None, force=False):

	obo = 'goslim_generic.obo' if 'slim' in ontology else 'go-basic.obo'

	namespace_filter = get_namespace_filter(exclude)

	if force & (os.path.isfile(obo)):
		os.remove(obo)

	obo_fname = goatools.base.download_go_basic_obo(obo)

	obodag = GODag(obo_fname)

	return {term_id:term.name for term_id,term in obodag.items() if namespace_filter(term.namespace)}

def get_gene_pathway_mapping(organism=9606, annotations=None, force=False):

	if force & (os.path.isfile('gene2go')):
		os.remove('gene2go')

	if annotations is None:
		annotation_filter = lambda x: True
	else:
		annotation_filter = lambda x: x in annotations

	gene2go = Gene2GoReader(goatools.base.download_ncbi_associations(), taxids=[organism]).get_id2gos_nss()

	gene_pathway_mapping = []

	for gene,gos in gene2go.items():
		gene_pathway_mapping.extend([(gene,go) for go in gos if annotation_filter(go)])

	return gene_pathway_mapping

def get_namespace_filter(exclude):

	namespace_map = {'BP':'biological_process','biological_process':'biological_process',
				 'MF':'molecular_function','molecular_function':'molecular_function',
				 'CC':'cellular_component','cellular_component':'cellular_component'
				}

	if exclude is None:
		return lambda x: True

	else:
		if isinstance(exclude,str):
			namespace = [exclude]

		try:
			namespace = [namespace_map[ns] for ns in exclude]
		except:
			raise ValueError(f'`exclude`={exclude} not recognized. Acceptable values are {namespace_map.keys()}')

		return lambda x: x not in namespace