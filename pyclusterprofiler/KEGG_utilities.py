import sharepathway
import re

__all__ = [
	'get_gene_id_mapping',
	'get_pathway_mapping',
	'get_gene_pathway_mapping'
	]

def get_gene_id_mapping(organism='hsa',force=False):
	return geneIDconv(organism=organism,force=force)

def get_pathway_mapping(organism='hsa', exclude=None, force=False):
	
	pathway_filter = get_pathway_filter(exclude,organism)

	pathway_tuples=sharepathway.parse_kegg.Request('list','pathway',organism, force=force)

	return {pathway_id:name for pathway_id,name in pathway_tuples if pathway_filter(pathway_id)}

def get_gene_pathway_mapping(organism='hsa', annotations=None, force=False):

	if annotations is None:
		annotations_filter = lambda x: True
	else:
		annotations_filter = lambda x: x in annotations

	return [(g,p) for g,p in sharepathway.parse_kegg.Request('link','path',organism, force=force) if annotations_filter(p)]

def get_pathway_filter(exclude=None,organism='hsa'):

	exclude_groups = ['human_diseases','organismal_systems']

	if exclude is None:
		return lambda x: True
	else:
		if organism != 'hsa':
			raise ValueError(f'Only exluding human (`organism`="hsa") KEGG pathways is implemented')

		if isinstance(exclude,str):
			exclude = [exclude]

		if not all([e in exclude_groups for e in exclude]):
			raise ValueError(f'`exclude`={exclude} not recognized. Acceptable values are {exclude_groups}')


		# this is hard-coded currently, since no clear way to separate KEGG pathway categories
		patterns = []
		if 'human_diseases' in exclude:
			patterns.extend([r'.*hsa05.*',r'.*hsa0493[01234].*',r'.*hsa049[45]0.*',r'.*hsa015.*'])
		if 'remove_organismal' in exclude:
			patterns.extend([r'.*hsa046[12457].*',r'.*hsa04062.*',r'.*hsa049[1267].*',r'.*hsa04935.*',
				r'.*hsa03320.*',r'.*hsa042[67].*',r'.*hsa047.*',r'.*hsa043[268].*',r'.*hsa0421[123].*'])

		return lambda x: not any([re.match(p,x) for p in patterns])

def geneIDconv(organism='hsa',force=False):
    # KEGG ID conv from NCBI-GeneID to KEGG ID , though most of them are the same
    keggd=sharepathway.parse_kegg.Request('conv',organism,'ncbi-geneid',force=force)
    Kgid = {}
    for line in keggd:
        Kgid[line[0].split(':')[1]] = line[1]
    return Kgid