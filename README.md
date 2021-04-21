# pyclusterprofiler

[![PyPI](https://img.shields.io/pypi/v/pyclusterprofiler.svg?color=green)](https://pypi.org/project/pyclusterprofiler)
[![Python Version](https://img.shields.io/pypi/pyversions/pyclusterprofiler.svg?color=green)](https://python.org)

A limited python implementation of [clusterProfiler] from R, borrowing some functions and concepts from [sharepathway] and [goatools].

Currently KEGG and GO interfaces are implemented.

----------------------------------

## Installation

You can install `pyclusterprofiler` via [pip]:

    pip install pyclusterprofiler

## Usage

	import pyclusterprofiler

To find enriched KEGG pathways in groupings ("cluster" column) of genes ("gene_id" column) identified in `df`:

	df_enrichment = pyclusterprofiler.compare_clusters(df,'cluster',database='KEGG')

Or using GO terms (instead using `database`="GO-slim" here will use reduced set of terms):
	
	df_enrichment = pyclusterprofiler.compare_clusters(df,'cluster',database='GO')

Example filter for any pathways/annotations with significant enrichment:
	
	significant_pathways = (df_enrichment
		.query('(corrected_pvalue<0.05)&(cluster_pathway_genes>3)')
		['pathway']
		.unique()
		)

Plot results as a dot plot:

	ax = pyclusterprofiler.dotplot(df_enrichment.query('pathway in @significant_pathways'))

### `compare_clusters` arguments

| argument | description |
|----------|-------------|
| `df` | dataframe with "gene_id" column containing NCBI gene id's and a column specifying group membership|
| `grouping` | column or list of columns in `df` to use for group membership |
| `enrichment_threshold` | threshold on ratio of observed/expected gene counts for test to include in results (default 1) |
| `correction` | method for correcting p-values for multiple hypothesis testing, used as argument to `statsmodels.stats.multitest.multipletests` (default "fdr_bh") |
| `organism` | organism databases to download. GO uses NCBI taxid; for KEGG see their [organism list]	(default is human databases for each) |
| `database` | "KEGG", "GO", or "GO-slim" (default "KEGG") |
| `exclude` | pathway/annotation groupings to exclude. For KEGG, can be "human_diseases", "organismal_systems," or a list of both (see [KEGG pathways]). For GO, can be "molecular_function","biological_process", "cellular_component", or a list of one or more (can also use abbreviations "MF","BP","CC" respectively) (default None) |
| `force` | force fresh download of databases, otherwise uses previously downloaded files if found in the current working directory (default False) |
| `verbose` | If True, prints provided NCBI gene id's that could not be found in the database (default True) |

## Contributing

Contributions are very welcome.

## License

Distributed under the terms of the [MIT] license,
"pyclusterprofiler" is free and open source software.

## Issues

If you encounter any problems, please [file an issue] along with a detailed description.

[MIT]: http://opensource.org/licenses/MIT
[file an issue]: https://github.com/lukebfunk/pyclusterprofiler/issues
[pip]: https://pypi.org/project/pip/
[clusterProfiler]: https://github.com/YuLab-SMU/clusterProfiler
[sharepathway]: https://github.com/GuipengLi/SharePathway
[goatools]: https://github.com/tanghaibao/goatools
[organism list]: https://www.genome.jp/kegg/catalog/org_list.html
[KEGG pathways]: https://www.genome.jp/kegg/pathway.html