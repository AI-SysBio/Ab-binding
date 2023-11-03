import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
df = pd.read_csv('sample_repertoire_raw.csv', sep = '\t')
df.dropna(inplace=True)

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
    return reverse_complement

N=100


AAseq =  df['pep'].to_numpy()
ll = np.array([len(seq) for seq in AAseq])
nwhere = np.where(ll>100)[0]

nk = np.random.randint(0, len(nwhere), size = N)

nkeep = nwhere[nk]

new_df = pd.DataFrame()
new_df['seqID'] = np.arange(N)
new_df['Vgene'] = df['V_sub'].to_numpy()[nkeep]
new_df['Vgene'] = np.array([V.split('+')[0] for V in new_df['Vgene']])
new_df['Jgene'] = df['J_sub'].to_numpy()[nkeep]
new_df['Jgene'] = np.array([J.split('+')[0] for J in new_df['Jgene']])
Seqs = df['seq'].to_numpy()[nkeep]
Seqs_r = [reverse_complement(seq) for seq in Seqs]
new_df['seq'] = Seqs_r

new_df.index = new_df['seqID']
new_df.drop('seqID', axis=1, inplace=True)

new_df.to_csv('sample_repertoire.csv', sep = '\t')


print(new_df)

#print(nkeep)