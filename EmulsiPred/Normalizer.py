import PepUtils as pu
import argparse
import numpy as np
import pandas as pd
import os

# python Normalizer.py -s /home/projects/vaccine/people/s132421/Arbejde/Data/NormalizationValues/randompep.fsa -o /home/projects/vaccine/people/s132421/Arbejde/Data/NormalizationValues
# python Normalizer.py -s /home/projects/vaccine/people/s132421/Arbejde/Data/NormalizationValues/randompeptide1000.fsa -o /home/projects/vaccine/people/s132421/Arbejde/Data/NormalizationValues

parser = argparse.ArgumentParser()
parser.add_argument('-s', dest='sequence', default='', help='Sequence to make normalization with in fasta')
parser.add_argument('-o',dest = 'out_dir', default ='', help='Path to directory to save output in')

# Define the parsed arguments
args = parser.parse_args()
sequence_file = args.sequence
out_dir = args.out_dir

def read_fasta_file(sequence_file):
    with open(sequence_file,'r') as data_file:

        data = ''
        acc = ''
        for line in data_file:
            if line[0] == '>':
                acc = line[1:].strip()
            elif line[0] != '>':
                data = data+ line.strip()
    return(acc,data)

# Take out the sequence name of the peptide + the full sequence
acc, rand_seq = read_fasta_file(sequence_file)

# Define the length of the sequences and which scheme is used to score them
alpha_score = pu.norm_alpha(rand_seq, w1=7, w2=31, hydro=1)
beta_score = pu.norm_beta(rand_seq, w1=7, w2=31, hydro=1)
gamma_score = pu.norm_gamma(rand_seq, w1=7, w2=31, hydro=1)

# Turn the results above into dataframes to make them easier to work with
a_df = pd.DataFrame(columns=['Length','Mean','Std'])
for length in alpha_score:
    a_df = pd.concat([a_df, pd.DataFrame([[str(length[0]), str(np.mean(length[1])), str(np.std(length[1]))]], index=['ALPHA'],columns=['Length', 'Mean', 'Std'])])

b_df = pd.DataFrame(columns=['Length','Mean','Std'])
for length in beta_score:
    b_df = pd.concat([b_df,pd.DataFrame([[str(length[0]), str(np.mean(length[1])), str(np.std(length[1]))]], index=['BETA'],columns=['Length', 'Mean', 'Std'])])

g_df = pd.DataFrame(columns=['Length','Cut','Mean','Std'])
for length in range(len(gamma_score)):
    for cut in gamma_score[length][1]:
        g_df = pd.concat([g_df,pd.DataFrame([[gamma_score[length][0],cut[0], str(np.mean(cut[1])), str(np.std(cut[1]))]], index=['GAMMA'],columns=['Length', 'Cut', 'Mean', 'Std'])])

# Save the dataframes
a_df.to_csv(out_dir + '/a_norm.csv', encoding='utf-8')
b_df.to_csv(out_dir + '/b_norm.csv', encoding='utf-8')
g_df.to_csv(out_dir + '/g_norm.csv', encoding='utf-8')