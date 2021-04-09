import seaborn as sns

__all__ = [
	'dotplot'
	]

def dotplot(df,x='cluster',y='pathway',size='enrichment',hue='-log10_corrected_pvalue',**kwargs):
	ax = sns.scatterplot(data=df,
	                x=x,
	                y=y,
	                size=size,
	                hue=hue,
	                size_norm=(1,10),
	                hue_norm=(2,10),
	                palette=sns.color_palette('viridis',as_cmap=True),
	                linewidth=0.5,
	                edgecolor='black'
	               )
	# ax.set_xscale('log')
	ax.set_xticks(list(range(df[x].max()+1)))
	ax.grid()
	ax.legend(loc=(1.1,0))