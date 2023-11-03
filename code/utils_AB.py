    
import os, sys
import numpy as np
import pandas as pd
import subprocess

from Bio.PDB import PDBParser
import Levenshtein
#install with pip install python-Levenshtein




def get_CDRs(AAseq_list, numbering = 'chothia', return_dict=False, fasta_filename = 'AAseqs.fasta'):
    #https://github.com/oxpig/ANARCI    
    #numbering = 'chothia', 'imgt', ...
    all_CDRs = dict()
    i = 0
    
    make_FASTA(AAseq_list, filename = fasta_filename)

    from Bio import SeqIO
    from tools.parapred.parapred.full_seq_processor import get_CDR_simple
    
    for record in SeqIO.parse(fasta_filename, "fasta"):
        CDRs = get_CDR_simple(str(record.seq).replace('X','').replace('*',''), scheme=numbering)
        #print(record.id, CDRs)
        try:
            all_CDRs[record.id] = [CDRs['CDR1'].replace('-',''),CDRs['CDR2'].replace('-',''),CDRs['CDR3'].replace('-','')]
        except:
            all_CDRs[record.id] = ['','','']
            
        i += 1
        if i%10000 == 0:
            print(' ',i)
        
        
    if return_dict:    
        return all_CDRs
      
    else:
        CDRs = np.array(list(all_CDRs.values()))
        
        return CDRs
    
    
    
def get_paratopes(AAseq_list, fasta_filename = 'AAseqs.fasta', output_file = 'paratopes.txt'):
    
    
    """
    Note: - Parapred only runs on Linux
          - parapred requires Anarci
    
    Note:
        https://github.com/eliberis/parapred
        Parapred is using an old version of tensorflow which is quite annoying.
        TO run it, you need to create a conda environment with a lower version of python (such as 3.6)
        Depending on the one you use, you should adapt the code below for the parapred cmw line
        
        you can install with:
        # conda create -n parapred python=3.6
        # conda activate parapred
        # conda install -c bioconda anarci
        # pip install Cython
        # python3 setup.py install
        
    Note that you can do all of these in the Windows sublinux system:
    On windows, you can launch unbuntu and launch the python scrit with the following commandes
    
    # Run Ubuntu as ADMINISTRATOR
    # conda activate parapred
    # cd /mnt/YOUR_PATH
    # python3 run_example.py
    """
    
    nseq = len(AAseq_list)
    AAseq = []
    seqIDs = []
    for si,seq in enumerate(AAseq_list):
        if '*' not in str(seq):
            try:
                AA = seq.replace('X','').replace('.','').replace('*','')
                AAseq.append(AA)
                seqIDs.append(str(si).zfill(6))
            except:
                pass
    make_FASTA(AAseq,fasta_filename, seqIDs = seqIDs) 

    #parapred_dir = "parapred-master"
    def launch_parapred(input_fasta, output_file):
        cwd = os.getcwd()
        filepath = os.path.join(cwd,input_fasta)
        filepath = '"%s"' % filepath
        outputpath = os.path.join(cwd,output_file)
    
        if os.name == 'nt':
            print("WARNING: PARAPRED can only be run on Linux, no new files will be created.")
        else:
            print('Running Parapred...')
            foutput = open(outputpath, "w")
            #process = subprocess.Popen("python -W ignore ./parapred-runner.py cdr SGYTFNSSDM", stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            cmd = "python ./tools/parapred/parapred-runner.py fasta {0}".format(filepath)
            print(cmd)
            process = subprocess.Popen(cmd, stdout=foutput, stderr=subprocess.PIPE, shell=True)#,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = process.communicate(timeout=400)
            print(out, err)
            print("  Parapred returned",process.poll())
            foutput.close() 

        # Read output from parapred
        paratope_threshold = 0.67 #asproposed in the Parapred paper
        AA_list = ['C','D','S','Q','K','I','P','T','F','N','G','H','L','R','W','A','V','E','Y','Y','M']
        result_dict = dict()
        paratope_dict = dict()
        try:
            f = open(outputpath,"r") 
        except:
            print('ERROR: No readable parapred output files was found')
            sys.exit()
        seqID = None
        CDR = -1
        for line in f:
            line = line.replace('\n','')
            if len(line) > 0:
                if ">" in line:
                    seqID = line.replace(">","").replace(" ","")
                    result_dict[seqID] = [[],[]]
                    paratope_dict[seqID] = []
                elif "{" in line:
                    #result_dict[seqID][0] = line
                    #print(line)
                    pass
                    
                elif "#" in line:
                    CDR = (CDR+1) % 3
                    #print("CDR",CDR+1)
                elif seqID is not None:
                    if line[0] in AA_list:
                        AA, value = line.split(" ")
                        result_dict[seqID][0].append(AA)
                        result_dict[seqID][1].append(value)
                        if float(value) > paratope_threshold:
                            paratope_dict[seqID].append(AA)
                        else:
                            paratope_dict[seqID].append("X")
                    pass
                    #print(line)
        f.close() 
        paratopes = ["".join(paratope_dict[seqID]) for seqID in paratope_dict]
        seqIDs = [seqID for seqID in paratope_dict]
        return np.array(paratopes), np.array(seqIDs)
    
    paratopes_, IDs = launch_parapred(fasta_filename, output_file)
    indices_  = np.array(IDs).astype(int)
    
    paratopes = []
    for i in range(nseq):
        if i in indices_:
            ii = np.where(indices_ == i)[0][0]
            paratopes.append(paratopes_[ii])
        else:
            paratopes.append('NA')
    
    return np.array(paratopes)



def paratope_distance(paratope1, paratope2):
    
    def count_residues(paratope):
        count = 0
        for i in range(len(paratope)):
            if paratope[i] != "X":
                count += 1
        return count
    
    
    if len(paratope1) != len(paratope2):
        return 1
    else:
        dist = 0
        nres = len(paratope1)
        for i in range(nres):
            if paratope1[i] == "X" or paratope2[i] == "X":
                pass
            else:
                if paratope1[i] != paratope2[i]:
                    dist += 1
        
        distance = dist/min(count_residues(paratope1),count_residues(paratope2))
        return distance
    



def make_FASTA(seq_list, filename, seqIDs=None):
    query_sequences = ''
    for i in range(len(seq_list)):
        if seqIDs is None:
            query_sequences += '>%s\n' % (str(i).zfill(6))
        else:
            query_sequences += '>%s\n' % (seqIDs[i])
        query_sequences += str(seq_list[i])
        query_sequences += '\n'    
    
    text_file = open(filename, "w")
    text_file.write(query_sequences)
    text_file.close()  
    
    
def Norm_Hamming_dist(string1,string2):
    """
    Returns the Hamming distance between string1 and string2
    If string1 and string2 have different length, the end of the longest string is cropped
    """
    if len(string1)==0 and len(string2)==0:
        return 1
    
    else:
        H = 0
        for i in range(min(len(string1),len(string2))):
            if string1[i] != string2[i]:
                H += 1
                
        return H/min(len(string1),len(string2))

def Norm_Levenshtein_dist(string1, string2):
    """
    [ref] Yujian, Li, and Liu Bo. "A normalized Levenshtein distance metric." 
    IEEE transactions on pattern analysis and machine intelligence 29.6 (2007): 1091-1095.
    """
    try:
        if len(string1)==0 and len(string2)==0:
            return 1
        else:
            Lev = Levenshtein.distance(string1, string2)
            norm_lev = 2*Lev/(len(string1)+len(string2)+Lev)
            
            return norm_lev
        
    except:
        return np.nan
    
    
def get_AAseqs(ntseq_list):
    AAseq = []
    for seq in ntseq_list:
        AA = most_likely_AAseq(seq)
        if AA[-1] == '*':
            AA = AA[:-1]
        AAseq.append(AA)
        
    return AAseq    
    
    
    
def most_likely_AAseq(sequence_nt):
    AAseq1 = translate_RNA(sequence_nt, lag=0)
    AAseq2 = translate_RNA(sequence_nt, lag=1)
    AAseq3 = translate_RNA(sequence_nt, lag=2)
    n1 = AAseq1.count('*')
    n2 = AAseq2.count('*')
    n3 = AAseq3.count('*')
    
    if n1 <= n2 and n1 <= n3:
        AAseq = AAseq1
    elif n2 <= n1 and n2 <= n3:
        AAseq = AAseq2
    else:
        AAseq = AAseq3
        
    if AAseq[-1] == 'X':
        AAseq = AAseq[:-1]
        
    return AAseq

def translate_RNA(sequence_nt, lag=0):
    """
    :param sequence: (str) a DNA sequence string
    :return: (str) a protein string from the forward reading frame 1
    """
    
    codontable = {'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
                  'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
                  'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
                  'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
                  'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
                  'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
                  'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
                  'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
                  'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
                  'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
                  'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
                  'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
                  'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
                  'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
                  'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
                  'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
                  '---': '-',
                  }
    
    seq = sequence_nt[lag:].upper()
    prot = []
    
    for n in range(0, len(seq), 3):
        if seq[n:n + 3] in codontable:
            residue = codontable[seq[n:n + 3]]
        else:
            residue = "X"
    
        prot.append(residue)
      
    return "".join(prot)



def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
    return reverse_complement




def fill_AA_seqs(nt_seqs, Vgenes): 
    
    
    """
    Fill the initial gap in the FWR1 region with the unmutated germline of the V gene (Vroot).
    Input:
       seq = IMGT gapped VDJ sequence
       Vroot = IMGT root V gene
    """
    
    from Bio import SeqIO

    def get_seq(database, ID):
        #find the sequence corresponding to the ID in the list 
        #(note that ID can just be a substring of string contained in the database)
    
        if isinstance(ID, float):
            return np.nan
        
        #remove the blank space at the end of the ID
        if ID[-1] == ' ':
            ID = ID[:-1]
        
        IDs = [database[j].id for j in range(len(database))]
        n=0
        index = -1
        for text in IDs:
            if ID in text:
                index = n
                break
            else:
                n += 1
        
        if index == -1:
            print("  Warning,",ID,"is not in found the database, return nan")
            return np.nan
        
        else:
            return ''.join(database[index].seq)
        
        
    def align_and_merge(seq,Vroot):
        """
        Fill the initial gap in the FWR1 region with the unmutated germline of the V gene (Vroot).
        VDJ should be longer than V0.
        Input:
           seq = IMGT gapped VDJ sequence
           Vroot = IMGT root V gene
        """
        
        seq_new = list(seq)
        ngap = 0
        for i in range(len(seq)):
            if seq[i] == '.':
                ngap += 1
            else:
                break
            
        for i in range(ngap):
            seq_new[i] = Vroot[i]
        
        return ''.join(seq_new)
    
    
    filled_ntseqs = []
    IMGT_V = list(SeqIO.parse("IGHV.fasta", "fasta"))
    for i in range(len(nt_seqs)):
        Vgene = Vgenes[i]
        seq = nt_seqs[i]
        Vroot = get_seq(IMGT_V, Vgene)
        full_seq = align_and_merge(seq,Vroot)
        filled_ntseqs.append(full_seq)




## All functions below are about abligity




def build_paratopes_pip(df, pdb_dir, H_only=False, pip_dir = None):
    """
    - Predict the paratopes residues with parapred
    - Generate the structures with RepertoireBuilder
    - Create pip interaction files from the generated structures 

    Note: - Parapred and Abligity only runs on Linux
          - RepertoireBuilder requires manual use of the website
    """
    
    nseq = len(df)
        
    if os.name == 'nt':
        print("Parapred can only be run on Linux")
    
    else:
        #2 - Build pip and ht files with parapred
        if pip_dir is None:
            pip_dir = pdb_dir
        else:
            pass
        
        for pdb in os.listdir(pdb_dir):
            if not pdb.endswith('.pdb'): continue
            pdb_name = pdb.split('.')[0]
            pip_from_parapred(pdb,pdb_dir,pip_dir)
            pipprefix = os.path.join(pip_dir,pdb_name)
            process = subprocess.Popen('./tools/Abligity/hash.o -i {0}.pip -o {0}_F.ht -b 1.0 -x 7b > /dev/null'.format(pipprefix),stdout=subprocess.PIPE, shell=True)
            process.wait()
            
            
        
        
        if H_only: #Create pip and ht files with HV only
            for pip in os.listdir(pip_dir):
                if not pip.endswith('.pip'): continue
                if "HV" in pip: continue 
                print(pip)
                pipprefix = os.path.join(pip_dir,pip).split('.')[0]
                pip_path = pipprefix+".pip"
                with open(pip_path, 'r') as file:
                    lines = file.readlines()
                file.close()
                                   
                pipprefix_HV = pipprefix + "_HV"
                pip_path_new = pipprefix_HV + ".pip"
                with open(pip_path_new, 'w') as new_file:
                    for line in lines:
                        if "H/" in line.strip("\n") or "pymolstring" in line.strip("\n"):
                            new_file.write(line)
                new_file.close()   
                
                process = subprocess.Popen('./tools/Abligity/hash.o -i {0}.pip -o {0}.ht -b 1.0 -x 7b > /dev/null'.format(pipprefix_HV),stdout=subprocess.PIPE, shell=True)
                process.wait()
                

    if not os.path.exists(pip_dir):
        os.makedirs(pip_dir)
    #3 - return paratopes and pip file list
    pip_filelist = [os.path.join(pip_dir,f) for f in os.listdir(pip_dir) if ("pip" in f and not "HV" in f)]
    paratope_HC, paratope_LC = get_pip_seq(pip_filelist)
    
    
    #4 - reorder the results and find failedcomputation for consistant indexes
    paratope_HC_ordered = []
    paratope_LC_ordered = []
    pip_list_ordered = []
    AbLigity_index = []
    for i in range(nseq):
        found = False
        for j in range(len(pip_filelist)):
            pip = pip_filelist[j]
            if str(i).zfill(6) == os.path.basename(pip).split('.')[0]:
                found = True
                
                paratope_HC_ordered.append(paratope_HC[j])
                paratope_LC_ordered.append(paratope_LC[j])
                pip_list_ordered.append(pip)
                AbLigity_index.append(j)
                
        if found == False:
            print("    pip", i,"was not found")
            paratope_HC_ordered.append("NA")
            paratope_LC_ordered.append("NA")
            pip_list_ordered.append("NA")
            AbLigity_index.append(-1)
    
    return pip_list_ordered, AbLigity_index




def get_pip_seq(pip_list):
    
    def sec_from_pip(filename):
        pip_table =  pd.read_csv(filename, delimiter = ',')
        residues = pip_table['pymolstring'].values
        sequenceH = []
        sequenceL = []
        for res in residues:
            if res.split('/')[0] == "H":
                sequenceH.append(res.split('/')[2])
            elif res.split('/')[0] == "L":
                sequenceL.append(res.split('/')[2])
            
        para_H = ''.join([str(elem) for elem in sequenceH])
        para_L = ''.join([str(elem) for elem in sequenceL])
    
        return para_H, para_L

    para_Hseq = []
    para_Lseq = []
    for pip in pip_list:
        para_H, para_L = sec_from_pip(pip)
        para_Hseq.append(para_H)
        para_Lseq.append(para_L)
        
    return np.array(para_Hseq), np.array(para_Lseq)




def compute_Abligity_similarity(df, pip_dir, H_only=False):
    
    
    nseq = len(df)
      
    pip_filelist = [os.path.join(pip_dir,f) for f in os.listdir(pip_dir) if ("pip" in f and not "HV" in f)]
    paratope_HC, paratope_LC = get_pip_seq(pip_filelist)
    
    
    paratope_HC_ordered = []
    paratope_LC_ordered = []
    pip_list_ordered = []
    AbLigity_index = []
    ii=0
    for i in range(nseq):
        found = False
        for j in range(len(pip_filelist)):
            pip = pip_filelist[j]

            if str(i).zfill(6) == os.path.basename(pip).split('.')[0]:
                found = True
                
                pip_list_ordered.append(pip)
                paratope_HC_ordered.append(paratope_HC[j])
                paratope_LC_ordered.append(paratope_LC[j])
                AbLigity_index.append(ii)
                ii+=1
                
        if found == False:
            print("    pip", i,"was not found")
            paratope_HC_ordered.append("NA")
            paratope_LC_ordered.append("NA")
            pip_list_ordered.append("NA")
            AbLigity_index.append(-1)
    
    if False:
        pass
    else:
        
        #1 - Computing similarity
        print("Similarity calculations")
        if not H_only:
            
            if os.name == 'nt':
                print("AbLigity can only be run on Linux")
            else:
                pip_list = '{0}/ablist.txt'.format(pip_dir) 
                process = subprocess.Popen('ls {0}/*_F.ht > {1}'.format(pip_dir,pip_list),stdout=subprocess.PIPE, shell=True)
                process.wait()
                process = subprocess.Popen('./similarity.o -l {0} -o {1}/sim.txt -c'.format(pip_list,pip_dir),stdout=subprocess.PIPE, shell=True)
                process.wait()
                
            sim_matrix = pd.read_csv('{0}/sim.txt'.format(pip_dir),header=None).fillna(0).values
            np.fill_diagonal(sim_matrix, 1)
            
            
            # Recovering the indexes from the initial dataset
            nseq = len(AbLigity_index)
            new_sim_matrix = np.zeros((nseq,nseq))
            for i in range(nseq):
                for j in range(i,nseq): 
                    ii = AbLigity_index[i]
                    jj = AbLigity_index[j]
                    if ii != -1 and jj != -1: #-1 means that the structural modelling failed for that sequence
                        new_sim_matrix[i,j] = sim_matrix[ii,jj]
                        new_sim_matrix[j,i] = new_sim_matrix[i,j]
                        
            return new_sim_matrix, paratope_HC_ordered, paratope_LC_ordered
            
        
        else: #Computing similarity for Heavy chain only
            if os.name == 'nt':
                print("AbLigity can only be run on Linux")
            else:
                pip_list_HV = '{0}/ablist_HV.txt'.format(pip_dir) 
                process = subprocess.Popen('ls {0}/*_HV.ht > {1}'.format(pip_dir,pip_list_HV),stdout=subprocess.PIPE, shell=True)
                process.wait()
                process = subprocess.Popen('./tools/Abligity/similarity.o -l {0} -o {1}/sim_H.txt -c'.format(pip_list_HV,pip_dir),stdout=subprocess.PIPE, shell=True)
                process.wait()
        
            sim_matrix_HV = pd.read_csv('{0}/sim_H.txt'.format(pip_dir),header=None).fillna(0).values
            np.fill_diagonal(sim_matrix_HV, 1)
        
            # Recovering the indexes from the initial dataset
            nseq = len(AbLigity_index)
            new_sim_matrix_HV = np.zeros((nseq,nseq))
            for i in range(nseq):
                for j in range(i,nseq): 
                    ii = AbLigity_index[i]
                    jj = AbLigity_index[j]
                    if ii != -1 and jj != -1: #-1 means that the structural modelling failed for that sequence
                        new_sim_matrix_HV[i,j] = sim_matrix_HV[ii,jj]
                        new_sim_matrix_HV[j,i] = new_sim_matrix_HV[i,j]
        
            return new_sim_matrix_HV, paratope_HC_ordered, paratope_LC_ordered
        
        
        
def pip_from_parapred(pdb,pdb_dir,pip_dir):
    
    res_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
                'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
                'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
                'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'TYS': 'Y', 
                'MET': 'M', 'MSE': 'M'}
    
    residuegroups = {'D':'F','E':'F','H':'E','K':'E','R':'E', 
                     'I':'A','L':'A','M':'C','P':'A','V':'A', 
                     'N':'G','Q':'G','C':'C','S':'D','T':'D', 
                     'A':'A','G':'A','F':'B','W':'B','Y':'B'}
    
    def abres2pym(res,abc):
        if abc[0] == abc[1]: abc = abc[0]+abc[1].lower()
        return '{0}/{1}{2}/{3}'.format('H' if res.parent.id==abc[0] else 'L',res.id[1],res.id[2] if res.id[2]!=' ' else '',res_dict[res.resname])

    
    overwrite = True
    pdbfile = os.path.abspath(os.path.join(pdb_dir,pdb))
    outf = os.path.join(os.path.abspath(pip_dir),pdbfile.split('/')[-1].split('.')[0]+'.pip')
    if not overwrite and os.path.exists(outf): return 0;
	
    print(pdbfile)
    #renumber_models(pdbfile, pdbfile, "chothia")
    process = subprocess.Popen("python tools/parapred/parapred-runner.py pdb '{0}'".format(pdbfile), shell=True)
    process.wait(timeout=400)
    #outs, errs = process.communicate(timeout=400)
    #print(outs, errs)
    print(pdbfile,process.poll())
	
	
    # Parapred > 67% in B-factor
    assert pip_dir!="", "Undefined PIP file output"
	
    if not os.path.exists(pdbfile): return 1
    s = PDBParser().get_structure(pdb,pdbfile)
    paratope_res = {abres2pym(res,'HL'):res for res in s.get_residues() if res.resname in res_dict and 'CA' in res and res['CA'].bfactor>=67}
    df_forcsv = pd.DataFrame([[res, residuegroups[res.split('/')[-1]]]+['{:.3f}'.format(c) for c in paratope_res[res]['CA'].coord] for res in paratope_res],columns=['pymolstring','residuetypes','x','y','z'])
    df_forcsv.to_csv(outf)
    return 0