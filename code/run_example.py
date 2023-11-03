import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 20})

import os, sys


from utils_AB import Norm_Levenshtein_dist, paratope_distance


def main():
    
    
    #preprocess_data() 
    #Check out the script at the bottom to see how to preprocess your data
    #You need the CDRs of all your sequences, as well as the predicted paratopes and the structures if you wanna use all metrics
   
    
    df = pd.read_csv('processed_data/sample_repertoire_processed.csv', sep='\t', index_col = 0)
    nseq = len(df)
    
    compute_abligity_similarity(df)
    
    sys.exit()
    #Activate this script if you want to compute the Abligity similarity. You can only do it after you got all the required pdb files.
    #Note that this is quite slow to compute
    
    seq = df['AAseq'].to_numpy()
    pdb_list = df['pdb_path'].to_numpy()
    CDR_seq = [[df['CDR1'][i],df['CDR2'][i],df['CDR3'][i]] for i in range(nseq)]
    paratope = df['paratope'].to_numpy()
    
    
    
    
    
    
    #1 - Build pairwise distance matrix of all sequences in our sample_repertoire:
    metrics = ["Lev_seq","CDR_dist", "paratype", "AbLigity", "RMSD", "Tm-score"] #RMSD and Tm-score are quite slow
    
    metrics = ["Lev_seq","CDR_dist", "paratype", "AbLigity"] 
    

    if "RMSD" in metrics or 'Tm-score' in metrics:
        import tmscoring #requires newever version of numpy, incompatible with scipy 
        #NumPy version >=1.16.5 and <1.23.0 is required for this version of SciPy (detected version 1.24.2)
        
    if "AbLigity" in metrics:
        Ab_ligity_sim = np.load('processed_data/AbLigity_sim.npy')
    
    recompute_distance_matrix = True
    if recompute_distance_matrix: 
        dist_matrix = np.ones((nseq,nseq,len(metrics)))
        for metric in metrics:
            print('  Computing distance matrice for metric %s...' % metric)
            for i in range(nseq):
                pdb_i = pdb_list[i]
                
                for k in range(len(metrics)):
                    dist_matrix[i,i,k] = 0 
                
                for j in range(i+1,nseq):
                    pdb_j = pdb_list[j]
                    if metric == 'Lev_seq':
                        dist_matrix[i,j,0] = Norm_Levenshtein_dist(seq[i],seq[j])
                        dist_matrix[j,i,0] = dist_matrix[i,j,0]
                    
                    elif metric == 'CDR_dist':
                        nCDRs = len(CDR_seq[0])
                        dist_matrix[i,j,1] = np.sum([Norm_Levenshtein_dist(CDR_seq[i][k],CDR_seq[j][k]) for k in range(nCDRs)])/nCDRs
                        if np.isnan(dist_matrix[i,j,1]):
                            dist_matrix[i,j,1] = 0
                        dist_matrix[j,i,1] = dist_matrix[i,j,1]
                        
                        
                    elif metric == 'paratype-original': #will give a score to 0 to sequences with different CDRs
                        dist_matrix[i,j,2] = paratope_distance(paratope[i], paratope[j])
                        dist_matrix[j,i,2] = dist_matrix[i,j,2]
                        
                    elif metric == 'paratype':
                        dist_matrix[i,j,2] = Norm_Levenshtein_dist(paratope[i].replace('X',''), paratope[j].replace('X',''))
                        dist_matrix[j,i,2] = dist_matrix[i,j,2]

                    
                    elif metric == 'AbLigity':
                        dist_matrix[i,j,3] = 1 - Ab_ligity_sim[i,j]
                        dist_matrix[j,i,3] = dist_matrix[i,j,3]
                    
                    elif metric == 'RMSD':
                        if os.path.isfile(pdb_i) and os.path.isfile(pdb_j):
                            alignment = tmscoring.TMscoring(pdb_i, pdb_j)
                            dist_matrix[i,j,4] = alignment.rmsd(**alignment.get_current_values())
                            if dist_matrix[i,j,4] == 0:
                                dist_matrix[i,j,4] = 1
                        else:
                            dist_matrix[i,j,4] = 0
                        dist_matrix[j,i,4] = dist_matrix[i,j,4]
                    
                    elif metric == 'Tm-score':
                        if os.path.isfile(pdb_i) and os.path.isfile(pdb_j):
                            alignment = tmscoring.TMscoring(pdb_i, pdb_j)
                            dist_matrix[i,j,5] = 1 - alignment.tmscore(**alignment.get_current_values())
                            if dist_matrix[i,j,5] == 0:
                                dist_matrix[i,j,5] = 1
                        else:
                            dist_matrix[i,j,5] = 0
                        dist_matrix[j,i,5] = dist_matrix[i,j,5]
                    
        if "RMSD" in metrics:   
            dist_matrix[:,:,4] = dist_matrix[:,:,4]/np.nanmax(dist_matrix[:,:,4]) # normalize the RMSD
            
        np.save('Dist_matrix.npy', dist_matrix)
        
    else:
        dist_matrix = np.load('Dist_matrix.npy')
        
        
        
        
    #plot and cl
    for k in range(dist_matrix.shape[2]):
   
        M_dist = dist_matrix[:,:,k]
        plt.figure(figsize=(10,10))
        plt.title(metrics[k])
        plt.imshow(1-M_dist, cmap='inferno', vmin=0, vmax=1, origin='lower')
        plt.colorbar()
        plt.show()    
    
    
    
    
def compute_abligity_similarity(df):
    
    from utils_AB import build_paratopes_pip, compute_Abligity_similarity
    
    build_paratopes_pip(df, pdb_dir = 'processed_data/PDB_structures/RepBuilder/Sample_seq', pip_dir = 'processed_data/PIP_files')
    sys.exit()
    Ab_ligity_sim, _, _ = compute_Abligity_similarity(df, pip_dir = 'processed_data\PIP_files')
    np.save('processed_data/AbLigity_sim.npy', Ab_ligity_sim)
    


    
def preprocess_data():
    
    from utils_AB import get_AAseqs, get_CDRs, get_paratopes, make_FASTA
    
    #1 - Get CDRs and paratopes from anarci and parapred
    df = pd.read_csv('raw_data/sample_repertoire.csv', sep='\t', index_col = 0)
    nseqs = len(df)
    seqs = df['seq'].to_numpy()
    AALseqs = np.array(get_AAseqs(seqs))
    df['AAseq'] = AALseqs
    

    print('Identifying CDRs')
    CDRs = get_CDRs(AALseqs, fasta_filename = 'processed_data/AAseqs.fasta')
    
    print('Identifying Paratopes')
    paratopeL = get_paratopes(AALseqs, fasta_filename = 'processed_data/AAseqs.fasta', output_file = 'processed_data/paratopes.txt')
    CDRL1s = CDRs[:,0]
    CDRL2s = CDRs[:,1]
    CDRL3s = CDRs[:,2]
    
    df['CDR1'] = CDRL1s
    df['CDR2'] = CDRL2s
    df['CDR3'] = CDRL3s
    df['paratopeL'] = paratopeL 
    print(df)
    
    
    
    #-> prepare the file to submite to repertoire builder
    heavy_chain_seq = 'QVQLQQPGAELVKPGASVKMSCKASGYTFTSYNMHWVKQTPGQGLEWIGAIYPGNGDTSYNQKFKGKATLTADKSSSTAYMQLSSLTSEDSAVYYCARGEAMDYWGQGTSVTVSS'
    AAHseqs = [heavy_chain_seq for i in range(nseqs)]
    fastafile_h = 'processed_data/PDB_structures/rep_fasta_h.fasta'
    fastafile_l = 'processed_data/PDB_structures/rep_fasta_l.fasta'
    make_FASTA(AALseqs, fastafile_l)
    make_FASTA(AAHseqs, fastafile_h)
    #-> Now you need to submit these two fasta files manually to https://sysimm.org/rep_builder/

    
    pdb_paths = ['processed_data/PDB_structures/RepBuilder/Sample_seq/%s.pdb' % str(i).zfill(6) for i in range(nseqs)]
    df['pdb_path'] = pdb_paths
    
    
    df.to_csv('processed_data/sample_repertoire_processed.csv', sep='\t')
    
    print(df)

    
    
    


if __name__ == "__main__":
    main()