import os,json,ast
import numpy as np
import pandas as pd
from multiprocessing.pool import Pool
from Bio.PDB import *
from scipy.stats import itemfreq

from plot_functions import generate_plots
from utils import *

def run_main(args):
	########### Begin generating PIPs ###########
	parsed_pips = []
	if args.mode == 0:
		if args.v: print("Unbound structures")
		# Check if modelling is needed
		if len(os.listdir(args.pdb_dir)) == 0:
			if args.v: print("Need to run ABodyBuilder before this step")
			return 1;
		
		pool = Pool(args.num_procs)
		done = pool.starmap(pip_from_parapred, [(args,pdb) for pdb in os.listdir(args.pdb_dir) if pdb.endswith('.pdb')]) 
		pool.close()
		for pdb in os.listdir(args.pdb_dir):
			if not pdb.endswith('.pdb'): continue
			pdb_chains = pdb.split('.')[0]
			pip_from_parapred(args,pdb)
			pipprefix = os.path.join(args.pip_dir,pdb_chains)
			process = subprocess.Popen('./hash.o -i {0}.pip -o {0}.ht -b 1.0 -x 7b > /dev/null'.format(pipprefix),stdout=subprocess.PIPE, shell=True)
			process.wait()
			parsed_pips.append(pdb_chains)
	elif args.mode == 1:
		if args.v: print("Bound structures")
		# Generate hash table
		if args.v: print("Generate Hash Table")
		pdblist = [pdb.replace('.pdb','') for pdb in os.listdir(args.pdb_dir) if pdb.endswith('.pdb')]
		parsed_pips = hash_list(args,pdblist) 

		
	elif args.mode == 2:
		if args.v: print("Bound structures from SAbDab")
		try:
			from sabdab_data import writefvag
		except ImportError as e:
			print(e)
			print("Check your installation of these packages")
			sys.exit(1)
		
		args.input = sorted(args.input)
		args.num_input = len(args.input)
	
		# Write structures
		if args.v: print("Writing Structures")
		pool = Pool(args.num_procs)
		done = pool.map(writefvag, [(x,args) for x in args.input]) 
		pool.close()	

		# Generate hash table
		if args.v: print("Generate Hash Table")
		pdblist = [pdb.replace('.pdb','') for pdb in os.listdir(args.pdb_dir) if pdb.endswith('.pdb')]
		parsed_pips = hash_list(args,pdblist) 

			
	########### End Generating PIPs ###########
	
	########## Begin Comparing PIPs ###########
	if args.v: print("Similarity calculations")
	args.pip_list = '{0}/ablist.txt'.format(args.pip_dir)
	args.simcsv = '{0}/sim.csv'.format(args.pip_dir)
	
	print(os.getcwd())
	process = subprocess.Popen('ls {0}/*.ht > {1}'.format(args.pip_dir,args.pip_list),stdout=subprocess.PIPE, shell=True)
	process.wait()
	process = subprocess.Popen('./similarity.o -l {0} -o {1}/sim.txt -c'.format(args.pip_list,args.pip_dir),stdout=subprocess.PIPE, shell=True)
	process.wait()
	output_ablist = os.path.join(args.pip_dir,"ablist.json")
	with open(output_ablist,'w') as f:
		json.dump(args.input,f)
	aborder = pd.read_csv(args.pip_list,header=None,names=['Ab'])['Ab'].apply(lambda x: os.path.basename(x).split('.')[0]).tolist()
	sim_matrix = pd.read_csv('{0}/sim.txt'.format(args.pip_dir),header=None).dropna(axis=1,how='all').fillna(0).values
	np.fill_diagonal(sim_matrix, 1)

	sim_matrix = pd.DataFrame(sim_matrix,index=aborder,columns=aborder)
	sim_matrix.to_csv(args.simcsv)
	if args.v: print("Finished calculating similarity. File saved to {0}".format(args.simcsv))
	########### End Comparing PIPs ###########

	############ Begin Plotting ############
	generate_plots(args)
	############  End Plotting  ############
	return 0;
	
	
if __name__ == "__main__":
	import sys, os, subprocess
	import argparse

	parser = argparse.ArgumentParser(prog="Ab-Ligity") 

	parser.add_argument( '--mode',		  type=int, help="Modes: [0]: Unbound; 1: Bound; 2: Protein complexes from SAbDab", dest="mode", default=0)
	parser.add_argument( '--input','-i',  type=str, nargs="*", help="List of PDB_ABChainID_AGChainID",dest="input")
	parser.add_argument( '--seqs',		  type=str, help="CSV of sequences with Heavy and Light chain sequences (optional: affinity/cluster information)", dest="seq_list")
	parser.add_argument( '--pdb_dir',	  type=str, help="PDB Output directory",  								  dest="pdb_dir")
	parser.add_argument( '--pip_dir',	  type=str, help="PIP Output directory",  								  dest="pip_dir")
	parser.add_argument( '--parapred',    type=str, help="Parapred directory",  								  dest="parapred_dir",default='parapred-master')
	parser.add_argument( '--scheme',      type=str, help="Numbering Scheme",  								  dest="scheme", 	 	default="imgt")
	parser.add_argument( '--definition',  type=str, help="CDR definition",  								  dest="definition",	default="north")
	parser.add_argument( '--num_proc',    type=int, help="Number of cores",									  dest="num_procs", 	default=6)
	parser.add_argument( '--anchor_len',  type=int, help="Number of anchor residues to align",		  dest="anchor_len", default=5)
	parser.add_argument( '--binding_site',type=str, help="Binding site selection [paratope]/epitope", dest="binding_site",default="paratope",choices=['paratope','epitope'])
	parser.add_argument( '--affinity_label',type=str, help="Affinity Label", 						  dest="affinity_label",default="kD(nM)")
	parser.add_argument( '--overwrite',                     help="Overwrite existing", dest="overwrite", default=False)
	parser.add_argument( '--verbose','-v',			help="Verbose",										 	  dest="v",           action='store_true')
	args = parser.parse_args()
	
	# Standardise	
	args.pdb_dir = makedir(args.pdb_dir)
	args.pip_dir = makedir(args.pip_dir)  
	args.image_dir = makedir(os.path.join(args.pip_dir,'images'))
	args.parapred_dir = os.path.abspath(args.parapred_dir)
	if args.seq_list: args.seq_list = os.path.abspath(args.seq_list)  
	
	run_main(args)
