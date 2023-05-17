
import os,json,itertools,pickle,ast,operator,re,time
from ABDB import database as db
from ABDB.AbPDB.Select import select_all
from Bio.PDB import *
from anarci import run_anarci

from ABDB.AB_Utils.region_definitions import Accept
a = Accept()
#a.set_regions(['cdrs'])

class NotDisordered(select_all):
	def accept_residue(self,residue):
		if atom.resname not in res_dict: return 0
		else:	return 1
	def accept_atom(self,atom):
		atom.set_altloc(' ')
		return 1

def write_pdb(out_file,residuelist,forcechain=''):
	# Write to pdb
	pdboutf = open(out_file,'wb')
	n = 1
	for eachresidue in residuelist:
		if forcechain != '': eachresidue.parent.id=forcechain
		aapdb, n = eachresidue._get_output_string(NotDisordered(), n)
		pdboutf.write(aapdb)
	pdboutf.write('TER\nEND\n')
	pdboutf.close()
	
def writefvag((pdb_chains,args)):
	pdb,abchains,agchains = pdb_chains.split('_')
	
	out_file = os.path.join(args.pdb_dir,pdb_chains+'.pdb')	
	if os.path.exists(out_file): return pdb_chains
	
	# Quality criteria
	p = db.fetch(pdb)
	if args.v:
		"""
		Criteria:
		o Any resolution
		o X-ray structure
		o Paired
		o Complex / has protein antigen
		o No missing residues in the CDR
		"""
		rules = [pdb_chains]+ [False]*4
		if 'RAY' in p.method: rules[1] = True
		sabdab_abchains = abchains if abchains[0].upper() != abchains[1].upper() else abchains[0].upper()*2
		ab = p[sabdab_abchains]
		if ab.VH != 'NA' and ab.VL !='NA': rules[2] = True
		if p.complex and ab.antigen!=None and len(ab.antigen.get_sequence()) > 50: rules[3] = True
		a.numbering_scheme=args.scheme
		a.defintion=args.definition
		a.set_regions(['cdrs'])
		if True not in [a.accept(pos,chain) for chain,poss in ab.get_missing_residues(scheme=args.scheme).items() for pos in poss]: rules[4] = True
		if False in rules: print("{0}\nStructure solved by X-ray crystallography:{1}\nPaired Fv:{2}\nHas protein antigen:{3}\nComplete CDR loops:{4}\n".format(*rules)) 	
		
	residues = p.get_structure(scheme=args.scheme,definition=args.definition).get_residues()

	residuelist = []
	# Fv + Ag structure
	for r in residues:
		if (r.get_full_id()[3] in abchains and r.region != '?') or r.get_full_id()[3] in agchains:
			residuelist.append(r)
	write_pdb(out_file,residuelist)
	return pdb_chains
"""
def get_anchor(_numberedAb,chain):
	a.numbering_scheme = args.scheme
	a.definition = args.definition
	
	anchor = []
	for cdr_ in ['1','2','3']:
		cdr = chain + cdr_
		a.set_regions(['cdr'+cdr])
		cdrloop = [pos for pos in _numberedAb if a.accept(pos[0],cdr[-2]) and pos[1] != '-']
		if cdrloop == []: continue
		numberedAb = [pos for pos in _numberedAb if pos[1] != '-']
		# find 5 residues of anchors
		startposi = sorted(numberedAb).index(sorted(cdrloop)[0])
		endposi = sorted(numberedAb).index(sorted(cdrloop)[len(cdrloop)-1])
		anchor += [n for n,_ in sorted(numberedAb)[startposi-args.anchor_len:startposi]+sorted(numberedAb)[endposi:endposi+args.anchor_len]]
	return anchor

def get_anchor_residues(pdb,pdbprefix,abchains):
	parser = PDBParser()
	s = parser.get_structure(pdb,'{}.pdb'.format(pdbprefix))
	# Translation
	atoms = []
	for c in s[0]:
		if c.id not in abchains: continue
		if c.id == abchains[0]: ct = 'H'
		else: ct = 'L'
		
		sabdab_abchains = abchains if abchains[0].upper() != abchains[1].upper() else abchains[0].upper()*2
		anchors = get_anchor(db.fetch(pdb)[sabdab_abchains].get_numbering(scheme=args.scheme)[ct].items(),ct)

		for r in c.get_residues(): 
			if r.id[1:] in anchors:
				if 'CA' not in r: print pdb, r.id,r.resname,'error'
				atoms.append(r['CA'])
	return s, atoms
"""
