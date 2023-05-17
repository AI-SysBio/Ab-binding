import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
import matplotlib
from matplotlib import colors
import seaborn as sns
import sys, os
import tmscoring #requires newever version of numpy, incompatible with scipy 
  #NumPy version >=1.16.5 and <1.23.0 is required for this version of SciPy (detected version 1.24.2)
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster, set_link_color_palette
from scipy.spatial.distance import squareform
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score
from sklearn import preprocessing

from libAB import Norm_Levenshtein_dist, Norm_Hamming_dist, paratope_distance, translate_dna
from libAB import get_paratopes, build_paratopes_pip, compute_Abligity_similarity
from libEnsembleClustering import cluster_ensembles

from scipy import stats


recompute = False



def compute_paratopes(recompute = False):
    
    
    data = pd.read_csv('sequences.csv', sep = '\t', index_col = 0)
    print(data)
    print(data.index)
    Ab_ligity_sim, paratope2H, paratope2L = compute_Abligity_similarity(data, HV_only=False)
    
    np.save('Processed/Parapred/AbLigity_sim.npy', Ab_ligity_sim)

    
    
    if recompute:
        data = pd.read_csv('sequences.csv', sep = '\t', index_col = 0)
        #data['seqH'] = [s.replace('.','') for s in data['seqH_gapped']]
        #data['seqL'] = [s.replace('.','') for s in data['seqL_gapped']]
        
        """
        print(data)
        paratopes_H, paratopes_L, Hi, Li = get_paratopes(data)
        
        print(paratopes_H)
        print(paratopes_L)
        
        data['paratopeH'] = paratopes_H
        data['paratopeL'] = paratopes_L
        print(data)
        data.to_csv('sequences_new.csv', index=False, sep='\t')
        sys.exit()
        """
        
        data = pd.read_csv('sequences_new.csv', sep = '\t', index_col = 0)
        
        #pip_list_ordered, AbLigity_index = build_paratopes_pip(data, HV_only=False)
        Ab_ligity_sim = compute_Abligity_similarity(data, HV_only=False)
        np.save('Processed/Parapred/AbLigity_sim.npy', Ab_ligity_sim)
    
    else:
        Ab_ligity_sim = np.load('Processed/Parapred/AbLigity_sim.npy')
        
    return Ab_ligity_sim, paratope2H, paratope2L


def main():
    
    #perform_preliminary_steps()
    Ab_ligity_sim, paratope2H, paratope2L = compute_paratopes(recompute = False)
    
    #Compute distance matrix for all metrics:
    seq_df = pd.read_csv('sequences_new.csv', sep = '\t')
    pdb_dir = "Processed/Structures/RB"
    
    #Verify if all pdb were computed
    drop_indexes = []
    pdb_list = []
    keep_index = []
    for i in range(len(seq_df)):
        pdb_i = '%s/%s.pdb' % (pdb_dir,str(i).zfill(6))
        if not os.path.isfile(pdb_i):
            print(pdb_i,"was not found")
            drop_indexes.append(i)
            pdb_list.append(pdb_i)
        else:
            pdb_list.append(pdb_i)
            keep_index.append(i) 
            
    #seq_df = seq_df.drop(drop_indexes)
    nseq = len(seq_df)
    seq_labels = seq_df['sample'].to_numpy()
    seq_df.index = np.arange(nseq)
    VgeneH = seq_df['VgeneH'].to_numpy()
    JgeneH = seq_df['JgeneH'].to_numpy()
    VgeneL = seq_df['VgeneL'].to_numpy()
    JgeneL = seq_df['JgeneL'].to_numpy()
    CDR_seq = [[seq_df['CDRH1'][i],seq_df['CDRH2'][i],seq_df['CDRH3'][i],seq_df['CDRL1'][i],seq_df['CDRL2'][i],seq_df['CDRL3'][i]] for i in range(nseq)]
    seqH = seq_df['seqH_gapped'].to_numpy()
    seqL = seq_df['seqL_gapped'].to_numpy()
    paratopeH = seq_df['paratopeH'].to_numpy()
    paratopeL = seq_df['paratopeL'].to_numpy()
    metrics = ["Hamming_seq","CDR_dist", "paratype", "paratype2", "AbLigity", "RMSD", "Tm-score", "Qscore"]#, "Qscore"]  
    metrics = ["Hamming_seq","CDR_dist", "paratype", "paratype2", "AbLigity", "RMSD", "Tm-score"]
    #metrics = ["Hamming_seq","CDR_dist", "paratype", "paratype2", "AbLigity"]
    #metrics = ["Hamming_seq","CDR_dist", "RMSD", "Tm-score"]#, "Qscore"]  
    
    #metrics = ["CDR_dist", "paratype"]#, "Qscore"]
    #metrics = ["Hamming_seq","CDR_dist"]
    dist_matrix = np.ones((nseq,nseq,8))


    
        
    # 1 - Sequences based metrics: 
    if recompute: 
        for metric in metrics:
            print('  Computing distance matrice for metric %s...' % metric)
            for i in range(nseq):
                pdb_i = pdb_list[i]
                
                for k in range(len(metrics)):
                    dist_matrix[i,i,k] = 0 
                
                for j in range(i+1,nseq):
                    pdb_j = pdb_list[j]
                    if metric == 'Hamming_seq':
                        dist_matrix[i,j,0] = (Norm_Hamming_dist(seqH[i],seqH[j]) + Norm_Hamming_dist(seqL[i],seqL[j]))/2
                        dist_matrix[j,i,0] = dist_matrix[i,j,0]
                    
                    elif metric == 'CDR_dist':
                        dist_matrix[i,j,1] = np.sum([Norm_Levenshtein_dist(CDR_seq[i][k],CDR_seq[j][k]) for k in range(6)])/6
                        if np.isnan(dist_matrix[i,j,1]):
                            dist_matrix[i,j,1] = 0
                        
                        dist_matrix[j,i,1] = dist_matrix[i,j,1]
                        
                    elif metric == 'paratype':
                        #dist_matrix[i,j,2] = (paratope_distance(paratopeH[i], paratopeH[j]) + paratope_distance(paratopeL[i], paratopeL[j]))/2
                        dist_matrix[i,j,2] = (Norm_Levenshtein_dist(paratopeH[i].replace('X',''), paratopeH[j].replace('X','')) + Norm_Levenshtein_dist(paratopeL[i].replace('X',''), paratopeL[j].replace('X','')))/2
                        dist_matrix[j,i,2] = dist_matrix[i,j,2]

                    elif metric == 'paratype2':
                        dist_matrix[i,j,3] = (Norm_Levenshtein_dist(paratope2H[i], paratope2H[j]) + Norm_Levenshtein_dist(paratope2L[i], paratope2L[j]))/2
                        dist_matrix[j,i,3] = dist_matrix[i,j,3]
                    
                    elif metric == 'AbLigity':
                        dist_matrix[i,j,4] = 1 - Ab_ligity_sim[i,j]
                        dist_matrix[j,i,4] = dist_matrix[i,j,4]
                    
                    elif metric == 'RMSD':
                        if os.path.isfile(pdb_i) and os.path.isfile(pdb_j):
                            alignment = tmscoring.TMscoring(pdb_i, pdb_j)
                            dist_matrix[i,j,5] = alignment.rmsd(**alignment.get_current_values())
                            if dist_matrix[i,j,5] == 0:
                                dist_matrix[i,j,5] = 1
                        else:
                            dist_matrix[i,j,5] = 0
                        dist_matrix[j,i,5] = dist_matrix[i,j,5]
                    
                    elif metric == 'Tm-score':
                        if os.path.isfile(pdb_i) and os.path.isfile(pdb_j):
                            alignment = tmscoring.TMscoring(pdb_i, pdb_j)
                            dist_matrix[i,j,6] = 1 - alignment.tmscore(**alignment.get_current_values())
                            if dist_matrix[i,j,6] == 0:
                                dist_matrix[i,j,6] = 1
                        else:
                            dist_matrix[i,j,6] = 0
                        dist_matrix[j,i,6] = dist_matrix[i,j,6]
                    
                    elif metric == 'Qscore':
                        dist_matrix[i,j,7] = 1
                        dist_matrix[j,i,7] = dist_matrix[i,j,7]
        if "RMSD" in metrics:   
            dist_matrix[:,:,5] = dist_matrix[:,:,5]/np.nanmax(dist_matrix[:,:,5])
            
        np.save('Dist_matrix.npy', dist_matrix)
        
    else:
        dist_matrix = np.load('Dist_matrix.npy')
        
        
    #dist_matrix =  dist_matrix[keep_index,:,:][:,keep_index,:]
    nseq = dist_matrix.shape[0]
    
    for k in range(dist_matrix.shape[2]):
        for i in range(dist_matrix.shape[0]):
            for j in range(i+1, dist_matrix.shape[0]):
                if dist_matrix[i,j,k] == 0:
                    dist_matrix[i,j,k] = 1
                    dist_matrix[j,i,k] = 1
                    
        
    
    #cluster with CDR:
    #M_dist = dist_matrix[:,:,3]
    #clusters = cluster_HAC(M_dist, d_threshold = 0.2, n_clusters = 30, plot = False)
    #nsorted = np.argsort(clusters)
    #dist_matrix = dist_matrix[nsorted,:,:][:,nsorted,:]
    
    print('%s Antobodies' % len(pdb_list))
    metrics = ["Hamming_seq","CDR_dist", "paratype", "paratype2", "AbLigity", "RMSD", "Tm-score"]
    
    
    cluster_list = []
    #for k in range(dist_matrix.shape[2]):
    for k in range(dist_matrix.shape[2]-1):
        
        print(metrics[k])
        print( '%.2f' % (1-np.mean(dist_matrix[:,:,k])), '%.2f' % (1-np.median(dist_matrix[:,:,k])))
        
        if k == 0:
            d_threshold = 0.15
        elif k == 1:
            d_threshold = 0.3  #optimal is 1 - 0.774 on SabDab
        elif k == 2:
            d_threshold = 0.4
        elif k == 3:
            d_threshold = 0.4
        elif k == 4:
            d_threshold = 0.9     #optimal is 1 - 0.192 on Sabdab
        elif k == 5:
            d_threshold = 0.2
        elif k == 6:
            d_threshold = 0.15
        
        else:
            d_threshold = 0.1
        
        
        M_dist = dist_matrix[:,:,k]
        clusters = cluster_HAC(M_dist, d_threshold = d_threshold, plot = False)
        #clusters = cluster_HAC(M_dist, n_clusters = 25, plot = False)
        nsorted = np.argsort(clusters)
        
        print('%s clusters found' % len(set(clusters)))
        

        
        #M_dist_sorted = M_dist
        M_dist_sorted = M_dist[nsorted,:][:,nsorted]
        nwhere = np.where(np.sum(1-M_dist_sorted, axis=1) ==1 )[0]
        
        plt.figure(figsize=(10,10))
        if k == 3 and False:
            plt.imshow(M_dist_sorted<0.95, cmap='inferno', origin='lower')
        else:
            plt.imshow(1-M_dist_sorted, cmap='inferno', vmin=0, vmax=1, origin='lower')
        plt.colorbar()
        plt.show()    
        
        cluster_list.append(clusters)

        #print(clusters[nsorted])
        print()
        print()
        print()
        
        
    nclsut = dist_matrix.shape[2]-1
    cluster_sim_matrix = np.ones((nclsut, nclsut))
    for ki in range(nclsut):
        for kj in range(ki+1,nclsut):
            similarity = adjusted_rand_score(cluster_list[ki], cluster_list[kj])
            cluster_sim_matrix[ki,kj] = similarity
            cluster_sim_matrix[kj,ki] = similarity
            
    plt.figure(figsize=(10,10))
    plt.imshow(cluster_sim_matrix, cmap='inferno', vmin=0, vmax=1, origin='lower')
    plt.xticks(np.arange(len(metrics)), metrics, rotation=90)
    plt.yticks(np.arange(len(metrics)), metrics)
    plt.colorbar()
    plt.show() 
    
    
    cluster_cor_matrix = np.ones((nclsut, nclsut))
    for ki in range(nclsut):
        for kj in range(ki+1,nclsut):
            similarity = stats.pearsonr(dist_matrix[:,:,ki].reshape(-1), dist_matrix[:,:,kj].reshape(-1))[0]
            cluster_cor_matrix[ki,kj] = similarity
            cluster_cor_matrix[kj,ki] = similarity
            
    plt.figure(figsize=(10,10))
    plt.imshow(cluster_cor_matrix, cmap='inferno', vmin=0, vmax=1, origin='lower')
    plt.xticks(np.arange(len(metrics)), metrics, rotation=90)
    plt.yticks(np.arange(len(metrics)), metrics)
    plt.colorbar()
    plt.show() 
        
    
    
    #find pairs that are consistently high:
    pair_dict = dict()
    pair_dict_value = dict()
    pair_dict_sim = dict()
    nab_ligity_pair = 0
    nparatype_pair = 0
    for seqi in range(nseq):
        for seqj in range(seqi+1,nseq):
            for ki in range(nclsut):
                sim_val = 1-dist_matrix[seqi,seqj,ki]
                if cluster_list[ki][seqi] == cluster_list[ki][seqj]:
                    pname = '(%s,%s)' % (seq_labels[seqi],seq_labels[seqj])
                    #pname = '(%s,%s)' % (seqi,seqj)
                    if pname in pair_dict_value:
                        pair_dict_value[pname] += 1
                        pair_dict[pname].append(metrics[ki])
                        pair_dict_sim[pname].append(float('%.2f' % sim_val))
                    else:
                        pair_dict_value[pname] = 1
                        pair_dict[pname] = [metrics[ki]]
                        pair_dict_sim[pname] = [float('%.2f' % sim_val)]
                        
                        
                if ki == 4 and sim_val > 0.1 and False:
                    nab_ligity_pair += 1
                    print(pname)
                    print('  SameCDR3 length: %s' % (int(len(CDR_seq[seqi][2]) == len(CDR_seq[seqj][2]))))
                    print('  CDRsim: %.2f' % (1-dist_matrix[seqi,seqj,1]))
                    print('  paratype: %.2f' % (1-dist_matrix[seqi,seqj,2]))
                    print('  paratype2: %.2f' % (1-dist_matrix[seqi,seqj,3]))
                    print('  AbLigity: %.2f' % (1-dist_matrix[seqi,seqj,4]))
                    print('  VHgene: %s - %s' % (VgeneH[seqi], VgeneH[seqj]))
                    print('  JHgene: %s - %s' % (VgeneL[seqi], VgeneL[seqj]))
                    print('  VLgene: %s - %s' % (JgeneH[seqi], JgeneH[seqj]))
                    print('  JLgene: %s - %s' % (JgeneL[seqi], JgeneL[seqj]))
                    print('  paratopeH: %s - %s' % (paratopeH[seqi], paratopeH[seqj]))
                    print('  paratopeL: %s - %s' % (paratopeL[seqi], paratopeL[seqj]))
                    print('  paratopeHX: %s - %s' % (paratopeH[seqi].replace('X',''), paratopeH[seqj].replace('X','')))
                    print('  paratopeLX: %s - %s' % (paratopeL[seqi].replace('X',''), paratopeL[seqj].replace('X','')))
                    print('  paratopeH2: %s - %s' % (paratope2H[seqi], paratope2H[seqj]))
                    print('  paratopeL2: %s - %s' % (paratope2L[seqi], paratope2L[seqj]))
                    print('  CDRH1: %s - %s' % (CDR_seq[seqi][0], CDR_seq[seqj][0]))
                    print('  CDRH2: %s - %s' % (CDR_seq[seqi][1], CDR_seq[seqj][1]))
                    print('  CDRH3: %s - %s' % (CDR_seq[seqi][2], CDR_seq[seqj][2]))
                    print('  CDRL1: %s - %s' % (CDR_seq[seqi][3], CDR_seq[seqj][3]))
                    print('  CDRL2: %s - %s' % (CDR_seq[seqi][4], CDR_seq[seqj][4]))
                    print('  CDRL3: %s - %s' % (CDR_seq[seqi][5], CDR_seq[seqj][5]))
                    print()
                    
                if ki == 2 and sim_val > 0.5 and False:
                    nparatype_pair += 1
                    print(pname, sim_val)
                    print('%.2f' % (1-dist_matrix[seqi,seqj,1]), '%.2f' % (1-dist_matrix[seqi,seqj,2]), '%.2f' % (1-dist_matrix[seqi,seqj,3]))
                    print()
    print (nab_ligity_pair, 'abligity pairs')
    print (nparatype_pair, 'paratype pairs')
                        
                        
    best_pairs, vals = keywithmaxval(pair_dict_value)
    pair_df = pd.DataFrame()
    for pi,pname in enumerate(best_pairs):
        si,sj = pname.replace('(','').replace(')','').split(',')
        seqi = np.where(seq_labels == si)[0][0]
        seqj = np.where(seq_labels == sj)[0][0]
        sameCDR = int(len(paratopeH[seqi]) == len(paratopeH[seqj]))
        CDRsim = 1-dist_matrix[seqi,seqj,1]
        Absim = 1-dist_matrix[seqi,seqj,4]
        parasim1 = 1-dist_matrix[seqi,seqj,2]
        parasim2 = 1-dist_matrix[seqi,seqj,3]
        parasim = (parasim1 + parasim2)/2
        struc_sim = 1-dist_matrix[seqi,seqj,6]
        
        
        pair_dict[pname] = []
        pair_dict_sim[pname] = []
        Ab1 = si
        Ab2 = sj
        
        Cl_cond = sameCDR
        CDR_cond = CDRsim > 0.7
        para_cond = parasim > 0.65
        AB_cond = Absim > 0.15
        struc_cond = struc_sim > 0.85
        
        
        if Cl_cond:
            pair_dict[pname].append('sameCDRl')
            pair_dict_sim[pname].append('%.2f' % sameCDR)
        
        if CDR_cond:
            pair_dict[pname].append('CDRsim')
            pair_dict_sim[pname].append('%.2f' % CDRsim)
            
        if para_cond:
            pair_dict[pname].append('parasim')
            pair_dict_sim[pname].append('%.2f' % parasim) 
            
        if AB_cond:
            pair_dict[pname].append('Absim')
            pair_dict_sim[pname].append('%.2f' % Absim) 
            
        if struc_cond:
            pair_dict[pname].append('strucsim')
            pair_dict_sim[pname].append('%.2f' % struc_sim) 
        
        if Cl_cond and struc_cond:
            cond = len(pair_dict[pname]) >= 3
        else:
            cond = len(pair_dict[pname]) >= 2
        
        if cond:
            
                print('', pname)
                print('  SameCDR3 length: %s' % (sameCDR))
                print('  CDRsim: %.2f' % (1-dist_matrix[seqi,seqj,1]))
                print('  paratype: %.2f' % (1-dist_matrix[seqi,seqj,2]))
                print('  paratype2: %.2f' % (1-dist_matrix[seqi,seqj,3]))
                print('  AbLigity: %.2f' % (1-dist_matrix[seqi,seqj,4]))
                print('  tm-score: %.2f' % (1-dist_matrix[seqi,seqj,6]))
                print('  VHgene: %s - %s' % (VgeneH[seqi], VgeneH[seqj]))
                print('  JHgene: %s - %s' % (VgeneL[seqi], VgeneL[seqj]))
                print('  VLgene: %s - %s' % (JgeneH[seqi], JgeneH[seqj]))
                print('  JLgene: %s - %s' % (JgeneL[seqi], JgeneL[seqj]))
                #print('  paratopeH: %s - %s' % (paratopeH[seqi], paratopeH[seqj]))
                #print('  paratopeL: %s - %s' % (paratopeL[seqi], paratopeL[seqj]))
                #print('  paratopeHX: %s - %s' % (paratopeH[seqi].replace('X',''), paratopeH[seqj].replace('X','')))
                #print('  paratopeLX: %s - %s' % (paratopeL[seqi].replace('X',''), paratopeL[seqj].replace('X','')))
                print('  paratopeH: %s - %s' % (paratope2H[seqi], paratope2H[seqj]))
                print('  paratopeL: %s - %s' % (paratope2L[seqi], paratope2L[seqj]))
                print('  CDRH1: %s - %s' % (CDR_seq[seqi][0], CDR_seq[seqj][0]))
                print('  CDRH2: %s - %s' % (CDR_seq[seqi][1], CDR_seq[seqj][1]))
                print('  CDRH3: %s - %s' % (CDR_seq[seqi][2], CDR_seq[seqj][2]))
                print('  CDRL1: %s - %s' % (CDR_seq[seqi][3], CDR_seq[seqj][3]))
                print('  CDRL2: %s - %s' % (CDR_seq[seqi][4], CDR_seq[seqj][4]))
                print('  CDRL3: %s - %s' % (CDR_seq[seqi][5], CDR_seq[seqj][5]))
                print()
                print()
                
                
                #pair_df = pair_df.append([[pname, sameCDR, pair_dict[pname],pair_dict_sim[pname]]])
                if sameCDR:
                    sameCDR_ = 'Yes'
                else:
                    sameCDR_ = 'No'
                list_all = [Ab1, Ab2, sameCDR_, '%.2f' % CDRsim, '%.2f' % parasim, '%.3f' % Absim, '%.2f' % struc_sim]
                list_all.append('%s - %s' % (VgeneH[seqi], VgeneH[seqj]))
                list_all.append('%s - %s' % (VgeneL[seqi], VgeneL[seqj]))
                list_all.append('%s - %s' % (JgeneH[seqi], JgeneH[seqj]))
                list_all.append('%s - %s' % (JgeneL[seqi], JgeneL[seqj]))
                list_all.append('%s - %s' % (paratope2H[seqi], paratope2H[seqj]))
                list_all.append('%s - %s' % (paratope2L[seqi], paratope2L[seqj]))
                list_all.append('%s - %s' % (CDR_seq[seqi][0], CDR_seq[seqj][0]))
                list_all.append('%s - %s' % (CDR_seq[seqi][1], CDR_seq[seqj][1]))
                list_all.append('%s - %s' % (CDR_seq[seqi][2], CDR_seq[seqj][2]))
                list_all.append('%s - %s' % (CDR_seq[seqi][3], CDR_seq[seqj][3]))
                list_all.append('%s - %s' % (CDR_seq[seqi][4], CDR_seq[seqj][4]))
                list_all.append('%s - %s' % (CDR_seq[seqi][5], CDR_seq[seqj][5]))
                pair_df = pair_df.append([list_all])
            
    #pair_df.columns = ['pair','SamCDR3length', 'methods','values']
    pair_df.columns = ['Ab1', 'Ab2', 'Same CDR length', 'CDRsim' ,'PARAsim','Abligity','TM-score','VHgenes','JHgenes','VLgenes','JLgenes','paratopeH','paratopeL','CDRH1','CDRH2','CDRH3','CDRL1','CDRL2','CDRL3']
    pair_df.to_csv('convergent_pairs_final.csv',index=False, sep='\t')
    pair_df.to_excel('convergent_pairs_final.xlsx',index=False)
    print(pair_df)
            
        
    
    
    
    
    
    """
    plt.figure(figsize=(10,10))
    plt.imshow(1-dist_matrix[:,:,1], cmap='inferno', vmin=0, vmax=1, origin='lower')
    plt.colorbar()
    plt.show()    
    

    plt.figure(figsize=(10,10))
    plt.imshow(1-dist_matrix[:,:,1], cmap='inferno', vmin=0, vmax=1, origin='lower')
    plt.colorbar()
    plt.show()    
    
    plt.figure(figsize=(10,10))
    plt.imshow(1-dist_matrix[:,:,2], cmap='inferno', vmin=0, vmax=1, origin='lower')
    plt.colorbar()
    plt.show()    
    
    plt.figure(figsize=(10,10))
    plt.imshow(1-dist_matrix[:,:,3], cmap='inferno', vmin=0, vmax=1, origin='lower')
    plt.colorbar()
    plt.show()  
    """

    sys.exit()      
    
    
    #cluster antibodies for each metrics
    sterotype_class = seq_df['stereotype_class'].to_numpy()
    disease = seq_df['disease'].to_numpy()
    nclsut_max = 38
    cluster_results = []
    n_stereotyped = np.where(sterotype_class >= 0)[0]
    cluster_true = sterotype_class[n_stereotyped]
    cluster_score = np.zeros((len(metrics),nclsut_max))
    for mi,metric in enumerate(metrics):
        print(metric)
        M_dist = dist_matrix[:,:,mi]
        
        for nclust in range(nclsut_max):
            if nclust > 1:
                clusters = cluster_HAC(M_dist, n_clusters = nclust, plot = False)
                cluster_predicted = clusters[n_stereotyped]
                cluster_score[mi,nclust] = adjusted_rand_score(cluster_predicted, cluster_true)
                cluster_score[mi,nclust] = adjusted_mutual_info_score(cluster_predicted, cluster_true)
        
        best_nclust = np.argsort(cluster_score[mi,:])[-1]
        clusters = cluster_HAC(M_dist, n_clusters = best_nclust, plot = False)
        print("Optimal number of clusters in %s" % best_nclust)
        print('Best Adjusted Rand index is %.2f' % cluster_score[mi,best_nclust])
        
        
        
        
        cluster_results.append(clusters)
        
        
        sorted_index = np.argsort(clusters)
        M_dist = M_dist[sorted_index,:][:,sorted_index]
        clusters = cluster_HAC(M_dist, n_clusters = best_nclust, plot = True)
    
        le = preprocessing.LabelEncoder()
        label_index = le.fit_transform(sterotype_class[sorted_index])
        cmap = colors.ListedColormap(["black","green", "blue", "gold", "red","darkviolet","gray"])
        bounds=[0,1,2,3,4,5,6]
        norm = colors.BoundaryNorm(bounds, cmap.N)
        plt.figure(figsize=(10,1))
        #plt.imshow(label_index.reshape(1,-1), aspect='auto', cmap='hsv')
        plt.imshow(label_index.reshape(1,-1), interpolation='nearest', origin='lower', aspect='auto',cmap=cmap, norm=norm)
        plt.axis('off')
        plt.show()    
        
        label_index = le.fit_transform(disease[sorted_index])
        cmap = colors.ListedColormap(["purple","orange"])
        #bounds=[0,1,2,3,4,5,6]
        #norm = colors.BoundaryNorm(bounds, cmap.N)
        plt.figure(figsize=(10,1))
        #plt.imshow(label_index.reshape(1,-1), aspect='auto', cmap='hsv')
        plt.imshow(label_index.reshape(1,-1), interpolation='nearest', origin='lower', aspect='auto',cmap=cmap)
        plt.axis('off')
        plt.show()  
        
        label_index = le.fit_transform(clusters)
        cmap = colors.ListedColormap(["gold", "cyan", "red","darkviolet","green", "blue","gray", "deepskyblue", "lawngreen","teal","maroon","chartreuse","crimson","magenta","tan","aquamarine","lightpink","lime","yellow","rosybrown"])
        #bounds=[0,1,2,3,4,5,6]
        #norm = colors.BoundaryNorm(bounds, cmap.N)
        plt.figure(figsize=(10,1))
        #plt.imshow(label_index.reshape(1,-1), aspect='auto', cmap='hsv')
        plt.imshow(label_index.reshape(1,-1), interpolation='nearest', origin='lower', aspect='auto',cmap=cmap)
        plt.axis('off')
        plt.show()  
        
        
        
        
    #compare each result with the adjusted rand index
    cluster_results = np.array(cluster_results)
    rand_index_matrix = np.ones((len(metrics),len(metrics)))
    for mi in range(len(metrics)):
        for mj in range(mi+1,len(metrics)):
            rand_index_matrix[mi,mj] = np.abs(adjusted_rand_score(cluster_results[mi], cluster_results[mj]))
            rand_index_matrix[mj,mi] = rand_index_matrix[mi,mj]
            
    plt.figure(figsize=(9.5,8))
    plt.imshow(rand_index_matrix, origin='lower', cmap=plt.cm.get_cmap('inferno'), vmin = 0, vmax=1)
    cbar = plt.colorbar()
    cbar.set_label('Adjusted rand index')
    plt.xticks(np.arange(len(metrics)),metrics,rotation=90)
    plt.yticks(np.arange(len(metrics)),metrics)
    plt.show()
    
    plt.figure(figsize = (9,4))
    for mi,metric in enumerate(metrics):
        plt.plot(np.arange(2,nclsut_max), cluster_score[mi,2:], label = metric, color=sns.color_palette()[mi+2], lw=2)
    plt.ylabel('Adjusted Rand index')
    plt.xlabel('Number of clusters')
    plt.legend(loc = 'upper right')
    plt.ylim(0,1)
    plt.show()
    
    
    #combine the results with Consensus clustering
    """
    Cluster-based Similarity Partitioning Algorithm (CSPA).
    Meta-CLustering Algorithm (MCLA)
    Hybrid Bipartite Graph Formulation (HBGF)
    Non negative Matrix Factorization (NMF) based consensus clustering
    """
    method = 'hbgf' #'cspa', 'mcla', 'hbgf', 'nmf'
    cluster_results = np.array([cluster_results[0],cluster_results[0],cluster_results[0],cluster_results[0]])
    clusters_ce = cluster_ensembles(cluster_results,solver = method, verbose = True, nclass = 3)
    
    """
    sorted_index = np.argsort(clusters_ce)
    M_dist = np.mean(dist_matrix,axis = 2)
    M_dist = M_dist[sorted_index,:][:,sorted_index]
    clusters = cluster_HAC(M_dist, n_clusters = best_nclust, plot = True)
    """
    
    le = preprocessing.LabelEncoder()
    label_index = le.fit_transform(sterotype_class)
    cmap = colors.ListedColormap(["black","green", "blue", "gold", "red","darkviolet","gray"])
    bounds=[0,1,2,3,4,5,6]
    norm = colors.BoundaryNorm(bounds, cmap.N)
    plt.figure(figsize=(10,1))
    #plt.imshow(label_index.reshape(1,-1), aspect='auto', cmap='hsv')
    plt.imshow(label_index.reshape(1,-1), interpolation='nearest', origin='lower', aspect='auto',cmap=cmap, norm=norm)
    plt.axis('off')
    plt.show()    
        
    label_index = le.fit_transform(disease)
    cmap = colors.ListedColormap(["purple","orange"])
    plt.figure(figsize=(10,1))
    plt.imshow(label_index.reshape(1,-1), interpolation='nearest', origin='lower', aspect='auto',cmap=cmap)
    plt.axis('off')
    plt.show()  
        
    label_index = le.fit_transform(clusters_ce)
    cmap = colors.ListedColormap(["gold", "cyan", "red","darkviolet","green", "blue","gray", "deepskyblue","lime", "lawngreen","teal","maroon","chartreuse","crimson","magenta","tan","aquamarine","lightpink","yellow","rosybrown"])
    plt.figure(figsize=(10,1))
    plt.imshow(label_index.reshape(1,-1), interpolation='nearest', origin='lower', aspect='auto',cmap=cmap)
    plt.axis('off')
    plt.show()
    
    cluster_predicted = clusters_ce[n_stereotyped]
    ce_rand_score = adjusted_rand_score(cluster_predicted, cluster_true)
    
    print("  %s Consensus clusters" % (np.max(clusters_ce)+1))
    print("  Adjusted rand Score of %.2f" % ce_rand_score)
    
    print(cluster_results)
    print(clusters_ce)
    

    

def perform_preliminary_steps():
    
    """
    Preliminary step:
     0 - Run the script once to make a FASTA file for HV and LV sequences
    
     1 - Run the alignement through IMGT HighV-quest (http://www.imgt.org/HighV-QUEST/home.action), download the result
         These are necessary to get all the CDR informations, and get where there is missing part of the sequences
         
     2 - After Its done, Run the HV and LV sequences through repertoire Builder to get the structures (find AAseqs_HC.csv and AAseqs_LC.csv files). 
         If applicable, you need to first fill the missing FR1 region with the unmutated germiline of Vgene, or the predictions will fail
         
     3 - predict paratopes and make pip files (AbLigity and parapred only runs on LINUX)
    """
    
    if not os.path.exists('Processed'):
        os.makedirs('Processed')
    
    
    # 1 - prepare the fasta file for IMGT submission
    data_file = "Data/Richard-Lymphoma-Abs.xlsx"
    df = pd.read_excel(data_file)
    Chain = []
    for Vgene in df['IGV']:
        if 'K' in Vgene or 'L' in Vgene:
            Chain.append('LC')
        else:
            Chain.append('HC')
    df['chain'] = Chain
    
    df_HC = df[df['chain'] == 'HC']
    df_LC = df[df['chain'] == 'LC'] 
    df_HC.reset_index(inplace=True)
    df_LC.reset_index(inplace=True)
    
    HC_seqs = df_HC["VH/VL Seq"].to_numpy()
    #make_FASTA(HC_seqs,"Processed/HC_seqs.fasta", seqIDs = list(df_HC['Sample']))
    make_FASTA(HC_seqs,"Processed/HC_seqs.fasta")
    LC_seqs = df_LC["VH/VL Seq"].to_numpy()
    #make_FASTA(LC_seqs,"Processed/LC_seqs.fasta", seqIDs = list(df_LC['Sample']))
    make_FASTA(LC_seqs,"Processed/LC_seqs.fasta")
    
    samples = df_HC['Sample']
    
    
    
    # 2 - make a dataframewith all the sequences
    data_file_HC = "Processed/IMGT/HC/2_IMGT-gapped-nt-sequences.txt"
    data_file_LC = "Processed/IMGT/LC/2_IMGT-gapped-nt-sequences.txt"
    
    data_H = pd.read_csv(data_file_HC, sep = '\t')
    seq_H = data_H["V-D-J-REGION"]
    Vgene_H = np.array([data_H['V-GENE and allele'][i].split(' ')[1].split('*')[0] if not isinstance(data_H['V-GENE and allele'][i], float) else "NA" for i in range(len(data_H))])
    Jgene_H = np.array([data_H['J-GENE and allele'][i].split(' ')[1].split('*')[0] if not isinstance(data_H['J-GENE and allele'][i], float) else "NA" for i in range(len(data_H))])
    CDR1_H = np.array([translate_dna(data_H['CDR1-IMGT'][i]).replace('.','') if not isinstance(data_H['CDR1-IMGT'][i], float) else "NA" for i in range(len(data_H))])
    CDR2_H = np.array([translate_dna(data_H['CDR2-IMGT'][i]).replace('.','') if not isinstance(data_H['CDR2-IMGT'][i], float) else "NA" for i in range(len(data_H))])
    CDR3_H = np.array([translate_dna(data_H['CDR3-IMGT'][i]).replace('.','') if not isinstance(data_H['CDR3-IMGT'][i], float) else "NA" for i in range(len(data_H))])
    AA_data = pd.read_csv("Processed/IMGT/HC/5_AA-sequences.txt", sep = '\t')
    AAseq_H = AA_data["V-D-J-REGION"]

    data_L = pd.read_csv(data_file_LC, sep = '\t')    
    seq_L = data_L["V-J-REGION"]
    Vgene_L = np.array([data_L['V-GENE and allele'][i].split(' ')[1].split('*')[0] if not isinstance(data_L['V-GENE and allele'][i], float) else "NA" for i in range(len(data_L))])
    Jgene_L = np.array([data_L['J-GENE and allele'][i].split(' ')[1].split('*')[0] if not isinstance(data_L['J-GENE and allele'][i], float) else "NA" for i in range(len(data_L))])
    CDR1_L = np.array([translate_dna(data_L['CDR1-IMGT'][i]).replace('.','') if not isinstance(data_L['CDR1-IMGT'][i], float) else "NA" for i in range(len(data_L))])
    CDR2_L = np.array([translate_dna(data_L['CDR2-IMGT'][i]).replace('.','') if not isinstance(data_L['CDR2-IMGT'][i], float) else "NA" for i in range(len(data_L))])
    CDR3_L = np.array([translate_dna(data_L['CDR3-IMGT'][i]).replace('.','') if not isinstance(data_L['CDR3-IMGT'][i], float) else "NA" for i in range(len(data_L))])
    AA_data = pd.read_csv("Processed/IMGT/LC/5_AA-sequences.txt", sep = '\t')
    AAseq_L = AA_data["V-J-REGION"]
    
    
    seq_df = pd.DataFrame()
    seq_df['sample'] = samples
    seq_df['CDRH1'] = CDR1_H
    seq_df['CDRH2'] = CDR2_H
    seq_df['CDRH3'] = CDR3_H
    seq_df['CDRL1'] = CDR1_L
    seq_df['CDRL2'] = CDR2_L
    seq_df['CDRL3'] = CDR3_L
    seq_df['VgeneH'] = Vgene_H
    seq_df['VgeneL'] = Vgene_L
    seq_df['JgeneH'] = Jgene_H
    seq_df['JgeneL'] = Jgene_L
    seq_df['seqH_gapped'] = seq_H
    seq_df['seqL_gapped'] = seq_L
    seq_df['AAseqH'] = AAseq_H
    seq_df['AAseqL'] = AAseq_L
    
    print(seq_df)
    
    #print(seq_df)
    
    seq_df.to_csv('sequences.csv', sep = '\t')
    #3 - Make files for repertoire Builder submission
    make_FASTA(AAseq_H,"Processed/AAseqs_HC.fasta")
    make_FASTA(AAseq_L,"Processed/AAseqs_LC.fasta")
    
    
    sys.exit()
 
    

    """
    sys.exit()
    

    data_file = "SMZL_heavy_light_chains_10.xlsx"
    data = pd.read_excel(data_file)
    print(data.columns)
    seqs = data["Sequence(nt)"].to_numpy()
    make_FASTA(seqs,"Processed/seqs.fasta")
    AAseqs = data["Sequence(AA)"].to_numpy()
    chain = data["CHAIN"].to_numpy()
    VH_index = np.where(chain == 'VH')[0]
    VL_index = np.where(chain == 'VL')[0]
    
    AAseqs_HC = AAseqs[VH_index]
    AAseqs_LC = AAseqs[VL_index]
    
    make_FASTA(AAseqs_HC,"Processed/AAseqs_HC.fasta")
    make_FASTA(AAseqs_LC,"Processed/AAseqs_LC.fasta")

    # 1 = Send sequence list to repertoire builder, Add a light chain to it to perform the modeling.
    #(this part has to be done manually)
    pdb_dir = "Processed/Structures/RB"
    """
    
    
    

def make_FASTA(seqs, filename, seqIDs = None):
    query_sequences = ''
    for i in range(len(seqs)):
        if seqIDs is None:
            query_sequences += '>%s\n' % (str(i).zfill(6))
        else:  
            query_sequences += '>%s\n' % (str(i).zfill(6) + (":%s" % seqIDs[i]))
        query_sequences += seqs[i]
        query_sequences += '\n'    
    
    text_file = open(filename, "w")
    text_file.write(query_sequences)
    text_file.close()    
    
    
    
    
    
def cluster_HAC(M_dist, d_threshold=None, n_clusters=2, plot = False):
    
    #tuto = https://joernhees.de/blog/2015/08/26/scipy-hierarchical-clustering-and-dendrogram-tutorial/
    
    """MI is a distance matrix n*n between each element
       d_treshold is the maximum distance of elements within cluster"""
    
    f = M_dist.shape[0]
    
    #Making the matrix fully symetric (if its not the case already)
    for fi in range(0,f):
        for fj in range(0,f):
            M_dist[fi,fj] = M_dist[fj,fi]
        
    #---------------------------------- Plot distance matrix ----------------------------------#
    if plot:
        plt.figure(figsize=(9.5,8))
        plt.imshow(M_dist, origin='lower', cmap=plt.cm.get_cmap('inferno_r'), vmin=0, vmax=1)
        plt.xlabel('Index')
        plt.ylabel('Index')
        cbar = plt.colorbar()
        cbar.set_label('Distance') #, rotation=270)
        #plt.clim(0,1)
        plt.show()
    
    #Setting diagonal to one (mandatory, important for the squareform)
    for fi in range(0,f): 
        M_dist[fi,fi] = 0
    M_dist_ss = squareform(M_dist) #the HAC algorithm need the squareform as an input
    
    
    d_result = []
    
    
    #--------------------------- HAC clustering ----------------------------#
    if plot:
        print ("\n")
        print ("  Calculating clusters with hierarchical clustering...") 
            
        
    Z = linkage(M_dist_ss, 'complete')
    if d_threshold is not None:
        clusters = fcluster(Z, d_threshold, criterion='distance')
    else:  
        clusters = fcluster(Z, n_clusters, criterion='maxclust')
   

    #plot dendrogram
    if plot:
        n_clust = int(np.max(clusters))
        d_result.append(n_clust)
        print ("    For a threshold distance of d =", d_threshold, "there is", n_clust, "clusters")
            
        cm = plt.get_cmap('tab20')
        colors = cm.colors
        set_link_color_palette([matplotlib.colors.rgb2hex(rgb[:3]) for rgb in colors])
        plt.figure(figsize=(10, 5))
        fancy_dendrogram(Z, p=840, truncate_mode='lastp', color_threshold=0, above_threshold_color='k', annotate_above=20, no_labels=True, max_d=d_threshold)
        plt.show()
    
    return clusters



def fancy_dendrogram(*args, **kwargs):
    #a build in fucntion to plot better dendrogram
    max_d = kwargs.pop('max_d', None)
    if max_d and 'color_threshold' not in kwargs:
        kwargs['color_threshold'] = max_d
    annotate_above = kwargs.pop('annotate_above', 0)

    ddata = dendrogram(*args, **kwargs)

    if not kwargs.get('no_plot', False):
        plt.title('Hierarchical Clustering Dendrogram')
        #plt.xlabel('sample index or (cluster size)')
        plt.ylabel('distance')
        for i, d, c in zip(ddata['icoord'], ddata['dcoord'], ddata['color_list']):
            x = 0.5 * sum(i[1:3])
            y = d[1]
            if y > annotate_above:
                plt.plot(x, y, 'o', c=c)
                plt.annotate("%.3g" % y, (x, y), xytext=(0, -5),
                             textcoords='offset points',
                             va='top', ha='center')
        if max_d:
            plt.axhline(y=max_d, c='k')
    return ddata
    
    
    
    
    
    
    
    
    
    
    
def create_seq_file_for_repBuilder(data):
    
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
    
    
    AAseqs = []
    IMGT_V = list(SeqIO.parse("IGHV.fasta", "fasta"))
    for i in range(len(data)):
        Vgene = data['Vgene'][i]
        seq = data['seq'][i]
        Vroot = get_seq(IMGT_V, Vgene)
        full_seq = align_and_merge(seq,Vroot)
        AAseqs.append(translate_dna(full_seq).replace('X',''))

    make_FASTA(AAseqs,"structure_analysis/Dominant_clone_sequences.fasta")  
        
        
        
    
    
def keywithmaxval(d,n_best=None):
     """ Return the n_best best keys having the highest value of a dictionary d"""  
     v=np.array(list(d.values()))
     k=np.array(list(d.keys()))
     
     if n_best is not None:
         best_indexes = v.argsort()[::-1][:n_best]
     else:
         best_indexes = v.argsort()[::-1]
         
     best_keys = k[best_indexes] 
     best_values = v[best_indexes]
         
     return best_keys,best_values
    
    
    
    
    
    
    
    
    
    

    
    
if __name__ == "__main__":
    plt.rcParams.update({'font.size': 20})
    main()