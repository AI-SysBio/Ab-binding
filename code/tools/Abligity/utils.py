# Ab-Ligity functions
import os,json,ast,subprocess
import numpy as np
import pandas as pd
from multiprocessing.pool import Pool
from Bio.PDB import *

from scipy.stats import itemfreq


res_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
				'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
				'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
				'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'TYS': 'Y', 
				'MET': 'M', 'MSE': 'M'}
residuegroups = { 
    'D':'F','E':'F','H':'E','K':'E','R':'E', 
    'I':'A','L':'A','M':'C','P':'A','V':'A', 
    'N':'G','Q':'G','C':'C','S':'D','T':'D', 
    'A':'A','G':'A','F':'B','W':'B','Y':'B' 
} 


def res2pym(atom):
	return '{0}/{1}{2}/{3}/{4}'.format(atom.parent.parent.id,atom.parent.id[1],atom.parent.id[2] if atom.parent.id[2]!=' ' else '',res_dict[atom.parent.resname],atom.id)

def make_pip_bio(args,pdb_chains):
	"""
	Create PIP files from input PDB
	Binding site definition: Biopython 4.5 \AA NeighborSearch
	
	Parameters
	----------
	args : class
		A class object that has the following attributes:
			binding_site : str
				Binding site of interest: ["paratope"]/"epitope"
				Default paratope
			pdb_dir : str
				PDB file directory
			pip_dir : str
				Output PIP file directory
	pdb_chains : str
		Four-letter PDB ID, antibody chain IDs and antigen chain IDs separated by "_"
	Returns
	-------
	None
	"""
	pdb,abchains,agchains = pdb_chains.split('_')
	if args.binding_site=="paratope":
		chains1 = abchains
		chains2 = agchains
	elif args.binding_site=="epitope":
		chains1 = agchains
		chains2 = abchains
	else:
		print("Unrecognised binding site selection {0}".format(args.binding_site))
		return;
		
	if abchains[0] == abchains[1]: abchains = abchains[0]+abchains[1].lower() # SCFv; do this after fetching structures
	outf = os.path.join(args.pip_dir,'{0}.pip'.format(pdb_chains))
	if not args.overwrite and os.path.exists(outf): return;

	pdbprefix = os.path.join(args.pdb_dir,pdb_chains)
	s = PDBParser(QUIET=True).get_structure(pdb_chains,pdbprefix+'.pdb')

	bs_res = []
	bs_res_d = {}

	ns = NeighborSearch(list(s.get_atoms()))
	
	for at1,at2 in ns.search_all(4.5,level='A'):
		at1c = at1.parent.parent.id 
		at2c = at2.parent.parent.id 
		if at1.parent.resname not in res_dict or at2.parent.resname not in res_dict: continue

		if (at1c in chains1 and at2c in chains2):
			resat = res2pym(at1)

			#if at1.name in ['CA','C','N','O']: continue
			res = '/'.join(resat.split('/')[:-1])
			if 'CA' not in at1.parent: continue
			bs_res.append(resat)
			if res not in bs_res_d: bs_res_d.update({res:at1.parent})

		elif (at1c in chains2 and at2c in chains1):
			resat = res2pym(at2)

			#if at2.name in ['CA','C','N','O']: continue
			res = '/'.join(resat.split('/')[:-1])
			if 'CA' not in at2.parent: continue
			bs_res.append(resat)
			if res not in bs_res_d: bs_res_d.update({res:at2.parent})
			
	# Rank paratope  residues 
	rank = itemfreq(['/'.join(n.split('/')[:-1]) for n in set(bs_res)])
	
	# Take the C-alpha coordiantes for pips
	df_forcsv = pd.DataFrame([[res, residuegroups[res.split('/')[-1]]]+['{:.3f}'.format(c) for c in bs_res_d[res]['CA'].coord] for res,_ in rank],columns=['pymolstring','residuetypes','x','y','z'])
	df_forcsv.to_csv(outf)
	return 0

def abres2pym(res,abc):
	if abc[0] == abc[1]: abc = abc[0]+abc[1].lower()
	return '{0}/{1}{2}/{3}'.format('H' if res.parent.id==abc[0] else 'L',res.id[1],res.id[2] if res.id[2]!=' ' else '',res_dict[res.resname])


def pip_from_parapred(args,pdb):
	pdbfile = os.path.join(args.pdb_dir,pdb)
	outf = os.path.join(os.path.abspath(args.pip_dir),pdbfile.split('/')[-1].split('.')[0]+'.pip')
	if not args.overwrite and os.path.exists(outf): return 0;
	
	#print(pdbfile)
	renumber_models(pdbfile, pdbfile, "chothia")
	# Parapred modify on the input Chothia-numbered PDB file. Binding probability replaces B-factor, in percentage (i.e. > 67% yields good precision according to the original paper)
	cwd = os.getcwd()
	os.chdir(args.parapred_dir)
	process = subprocess.Popen("python3.6 -m parapred pdb {0}".format(pdbfile), shell=True)#,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	process.wait(timeout=60)
	#outs, errs = process.communicate(timeout=400)
	print(pdbfile,process.poll())
	
	os.chdir(cwd)
	
	# Parapred > 67% in B-factor
	assert args.pip_dir!="", "Undefined PIP file output"
	
	if not os.path.exists(pdbfile): return 1
	s = PDBParser().get_structure(pdb,pdbfile)

	paratope_res = {abres2pym(res,'HL'):res for res in s.get_residues() if res.resname in res_dict and 'CA' in res and res['CA'].bfactor>=67}
	
	df_forcsv = pd.DataFrame([[res, residuegroups[res.split('/')[-1]]]+['{:.3f}'.format(c) for c in paratope_res[res]['CA'].coord] for res in paratope_res],columns=['pymolstring','residuetypes','x','y','z'])
	
	df_forcsv.to_csv(outf)
	return 0
def makedir(raw_dir):
	abs_dir = os.path.abspath(raw_dir)
	if not os.path.exists(abs_dir): os.mkdir(abs_dir)
	return abs_dir

def ligityhash(pipprefix,binding_site):
	if binding_site == 'paratope':
		process = subprocess.Popen('./hash.o -i {0}.pip -o {0}.ht -b 1.0 -x 7b > /dev/null'.format(pipprefix),stdout=subprocess.PIPE, shell=True)
	elif binding_site == 'epitope':
		process = subprocess.Popen('./hash.o -i {0}.pip -o {0}.ht -b 1.0 -x 7b > /dev/null'.format(pipprefix),stdout=subprocess.PIPE, shell=True)
	else:
		print("Shouldn't reach here")
	process.wait()
		
def hash_list(args,pdblist):
	parsed_pips = []
	for pdb_chains in pdblist:
	
		return_code = make_pip_bio(args,pdb_chains)
		if return_code == 1: 
			print("Errors in making PIPs for {0}; removing from similarity calculation".format(pdb_chains))
			continue
		pipprefix = os.path.join(args.pip_dir,pdb_chains)
		ligityhash(pipprefix,args.binding_site)
		parsed_pips.append(pdb_chains)
	return parsed_pips
	
def renumber_models(input_pdb, output_pdb, scheme):
	"""
	Change the numbering to Chothia and chain IDs to H and L for Parapred
	"""
	from AbPDB.AntibodyParser import AntibodyParser
	from AbPDB.Select import select_all
	p = AntibodyParser(QUIET=True)
	p.set_numbering_method("anarci")
	p.set_numbering_scheme(scheme)
	s = p.get_antibody_structure( "model", input_pdb )
	with open( output_pdb , 'w' ) as outfile:
		outfile.write( s._get_output_string( select_all(), 1 )[0] )

	
	
	
	
