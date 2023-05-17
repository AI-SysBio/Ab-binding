import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import json, itertools, os

abmparathres,inmparathres = 0.2, 0.7
colours_cycle = ["#77aadd","#EE8866","#44BB99","#EEDD88","#DDDDDD","#FFAABB","#BBCC33"]

def create_clustermap(simmat,colours,figtitle,savefig=None):
	"""
	Create cluster map
	Parameters
	----------
	simmat : pandas.DataFrame
		Symmetrical similarity matrix with index and columns as the labels
		The diagonal should be 1
	colours : list
		Colours of the labels
	figtitle : str
		Figure title
	savefig : figure file name
		File name for the figure if desired
		Default: None (Not saving the figure)
	Returns
	-------
	None
	"""
	fig = sns.clustermap(simmat,
	                     row_colors=colours,
	                     col_colors=colours,
	                     vmin=0,vmax=1,
	                     figsize=(24,24)
	                    )
	fig.fig.suptitle(figtitle)
	plt.setp(fig.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
	plt.setp(fig.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
	if savefig: fig.savefig(savefig,dpi=100)
	return 0;


def generate_plots(args):
	# Ab-Ligity similarity
	abligity_sim = pd.read_csv(args.simcsv,index_col=0)	
	
	with open(args.pip_list) as f:
		piplist = [os.path.basename(x).split('.')[0] for x in f.readlines()]
	
	dtype = None
	if args.seq_list: 
		seq_df = pd.read_csv(args.seq_list)
		if args.affinity_label:
			affinity = seq_df[['Ab',args.affinity_label]].set_index('Ab').dropna()
			dtype = "numeric" if affinity[args.affinity_label].dtype==np.float64 else "categorical"
			if dtype == "categorical": colour_dict = dict(zip(affinity[args.affinity_label].drop_duplicates().tolist(),colours_cycle))
		else:
			affinity = seq_df
		ablist = list(set(piplist) & set(affinity.index))
		if args.affinity_label: affinity = affinity.loc[ablist].sort_values(args.affinity_label)
		ablist = affinity.index.tolist() # do it again to get the sorted indices
		abligity_sim = abligity_sim[ablist].reindex(ablist)
	else:
		ablist = set(piplist)

	# Colours
	if dtype == "categorical":
		colours=[colour_dict[cluster] for ab in ablist for cluster in colour_dict if affinity[args.affinity_label][ab] == cluster]

	else:
		colours = ["tab:pink"]*len(ablist)
	
	
	########### Cluster / Heatmap plots	########### 

	# Raw read
	# Ab-Ligity
	create_clustermap(abligity_sim,
					colours,
					figtitle='Ab-Ligity',
					savefig=os.path.join(args.image_dir,'ab_cl_raw.png'))
	
	#Thresholding
	# Ab-Ligity
	create_clustermap(abligity_sim>=abmparathres,
					colours,
					figtitle='Ab-Ligity$\geq${0:.2f}'.format(abmparathres),
					savefig=os.path.join(args.image_dir,'ab_cl_threshold.png'))


	# Affinity and Structural similarity
	if not args.affinity_label: return 0


	###### Network plots ########

	# Continuous
	if dtype == "numeric":
		
		G = nx.from_numpy_matrix(abligity_sim.values,create_using=nx.Graph())
		G = nx.relabel_nodes(G,dict(enumerate(ablist)))
		
		
		fig,ax = plt.subplots(figsize=(12,12))
		pos = nx.spring_layout(G,seed=1)
		nx.draw_networkx_edges(G, pos, width=0.1, alpha=0.5)
		nx.draw_networkx_nodes(G,pos,
                       node_color=affinity[args.affinity_label],
                       with_labels=True,
                       labels=ablist,
                       cmap=plt.cm.Blues_r,
                       alpha=0.5
                       )
		nx.draw_networkx_labels(G,pos,labels=dict(zip(ablist,zip(ablist,affinity[args.affinity_label]))))
		
		
		
		#fig,ax = plt.subplots(figsize=(12,12))
		#nx.draw(G, with_labels=True,node_color=affinity[args.affinity_label], cmap=plt.cm.Blues)
		fig.savefig(os.path.join(args.image_dir,'ab_network_cont.png'),dpi=100)
	elif dtype == "categorical":
		# Categorical
		
		
		G = nx.from_numpy_matrix(abligity_sim.values,create_using=nx.Graph())
		G = nx.relabel_nodes(G,dict(enumerate(ablist)))
		
		fig,ax = plt.subplots(figsize=(12,12))
		pos = nx.spring_layout(G,seed=1)
		
		# plot edge
		nx.draw_networkx_edges(G, pos, width=0.1, alpha=0.5)
		# plot each category of nodes
		for cat,colour in colour_dict.items():
		    nodelist = affinity.index[affinity[args.affinity_label]==cat].tolist()
		    nx.draw_networkx_nodes(G,pos,
		                           nodelist=nodelist,
		                           node_color=colour_dict[cat],
		                           with_labels=True,
		                           label=cat
		                          )
		    nx.draw_networkx_labels(G,pos,labels=dict(zip(nodelist,nodelist)))
		plt.legend()
		fig.savefig(os.path.join(args.image_dir,'ab_network_cat.png'),dpi=100)

	return 0




