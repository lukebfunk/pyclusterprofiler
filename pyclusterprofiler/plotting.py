import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from copy import copy

__all__ = [
	'dotplot'
	]

def dotplot(df,x='cluster',y='pathway',size='log2_observed/expected',hue='-log10_corrected_pvalue',
			figsize=plt.rcParams['figure.figsize'],size_norm=(1,10),hue_norm=(1,10),**kwargs):

	df = df.astype({x:str})

	if 'ax' in kwargs.keys():
		ax = kwargs.pop('ax')
	else:
		fig, ax = plt.subplots(1,1,figsize=figsize)

	sns.scatterplot(data=df,
	                x=x,
	                y=y,
	                size=size,
	                hue=hue,
	                size_norm=size_norm,
	                hue_norm=hue_norm,
	                palette=sns.color_palette('viridis',as_cmap=True),
	                linewidth=0.5,
	                edgecolor='black',
	                ax=ax,
	                **kwargs
	               )

	hdl,_ = ax.get_legend_handles_labels()

	legend_colors = sns.color_palette('viridis',as_cmap=True)(np.linspace(0,255,5,dtype=int))
	legend_color_vals = np.linspace(*hue_norm,5)
	legend_color_header=hdl[0]

	if 'sizes' in kwargs.keys():
		sizes = kwargs['sizes']
		if isinstance(sizes,dict):
			sizes = sizes.values()
	else:
		sizes = np.r_[.5, 2] * np.square(plt.rcParams["lines.markersize"])
	
	min_max = min(sizes),max(sizes)

	legend_sizes = np.linspace(*min_max,5)
	legend_size_vals = np.linspace(*size_norm,5)
	legend_size_header = copy(hdl[0])
	legend_size_header.set_label(size)

	legend_elements = ([legend_color_header]+[plt.scatter([],[],marker='o',s=legend_sizes[3],color=c,
		                                     linewidth=0.5,edgecolor='k',label=str(cl)) 
	                              for c,cl in zip(legend_colors,legend_color_vals)]
	                   +[legend_size_header]+[plt.scatter([],[],marker='o',s=s,color='k',
	                   					label=str(sl)) 
	                              for s,sl in zip(legend_sizes,legend_size_vals)]
	                  )

	y_units = ax.yaxis.get_units()._mapping
	x_units = ax.xaxis.get_units()._mapping

	if max(y_units.values())<10:
		y_offset = max(y_units.values())*0.05
	else:
		y_offset = 0.5	

	for _,s in df.iterrows():
		t = ax.text(s[x],y_units[s[y]]-y_offset,int(s['cluster_pathway_genes']),ha='center',va='center')

	for _,s in df.drop_duplicates(x).iterrows():
		ax.text(s[x],ax.get_ylim()[1]-(y_offset),int(s['cluster_pathway_genes']/s['cluster_pathway_ratio']),
			ha='center',va='center',rotation=-90
			)

	ax.grid()
	ax.tick_params(axis='x',rotation=-90)
	ax.legend(handles=legend_elements,loc=(1.05,0.25))
	ax.text(np.mean(ax.get_xlim()),ax.get_ylim()[1]-(y_offset*2.5),'# of genes in cluster',ha='center')

	return ax