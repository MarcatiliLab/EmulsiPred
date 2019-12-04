from __future__ import print_function
import PepUtils as pu
import argparse
import operator
import numpy as np
import os
import pandas as pd

# python Gamma_emulsifier.py -s /home/projects/vaccine/people/s132421/Arbejde/Data/PeptideSequences.fsa -v /home/projects/vaccine/people/s132421/Arbejde/Data/NormalizationValues/g_norm.csv -o /home/projects/vaccine/people/s132421/Arbejde/Data/Results

# Parse it in
parser = argparse.ArgumentParser()
parser.add_argument('-s', dest='sequences', default='', help='Fasta file with all the sequences')
parser.add_argument('-v', dest='norm_vals', default='', help='csv file with mean and standard deviation values')
parser.add_argument('-o', dest='out_dir', default='', help='Output directory path')

# Define the parsed arguments
args = parser.parse_args()
sequences_fasta = args.sequences
norm_vec_file = args.norm_vals
out_dir = args.out_dir


def g_emul(sequences,norm_df):

    Kaa = pu.hydrophobicity_scale(1)
    g_best = []

    for name in sequences:
        g_scores={}
        seq = sequences[name]

        for i,res in enumerate(seq):

            for j in range(7, 31):
                if i + j <= len(seq):

                    # Calculate the hydrophobicity for the peptide
                    gamma = pu.g_emul_calcer(seq[i:(i+j)], Kaa)
                    gamma = sorted(gamma, key=operator.itemgetter(0), reverse=True)[0:10]

                    # Pick the right values for normalization
                    temp_df = norm_df.loc[(norm_df['Length'] == len(seq[i:(i + j)])) & (norm_df['Cut'] == gamma[0][1])]

                    # Pick cut which gives the highest score (used to represent that peptide) and normalize it
                    g_scores[name, seq[i:(i + j)], gamma[0][1]] = pu.z_normalize(gamma[0][0], temp_df.iloc[0][2], temp_df.iloc[0][3])

        # Make sure it is always only the 10000 best which are saved (makes it run faster)
        g_best = pu.highest_in_list(g_best, g_scores, 10000)

    # Groups repetitive sequences
    switched_gamma = pu.switch_dictionary_gamma(g_best)

    return(switched_gamma)

def cut_offs(results):
    # Removing peptides with less than 3 as a z-score and seen in less than 4 accession numbers
    cut_results = {}
    for i in results:
        if len(results[i][1]) >= 4 and results[i][0] >= 3:
            cut_results[i] = results[i]

    return (cut_results)

def counter(gamma):
    # Counts the total charge of each peptide
    total = []
    for i in gamma:
        neg = 0
        pos = 0
        for j in range(len(i[0])):
            if i[0][j] == 'D' or i[0][j] =='E':
                neg += 1
            elif i[0][j] == 'R' or i[0][j] == 'K':
                pos += 1

        diff = pos - neg
        total.append(diff)

    return(total)

def txt_file(results, charge, out_dir):
    # Writes the results in a nice viewable file
    with open(out_dir+'/g_results.txt', 'w') as outfile:

        outfile.write('-----' + '\t' + 'Charge' + '\t' + 'Cut' '\t' + 'Kyte-Doolittle' + '\t' + 'Sequence' + '\t\t\t' + 'Diff seqs' + '\n')
        counter = 0
        for i in results:
            outfile.write("GAMMA" + '\t' + str(charge[counter])  + '\t' + str(i[1]) +
                          '\t'+ '{:.5f}'.format(results[i][0]) + '\t\t' + "{:<31}".format(i[0]) + '\t' + str(len(results[i][1])) + '\t' + str(results[i][1]) + '\n')

            counter += 1

def fasta_file(results, charge, out_dir):
    # Fasta file used for clustering
    with open(out_dir+'/g_results.fsa', 'w') as outfile:
        counter = 0
        for i in results:
            outfile.write(">GAMMA|" + str(charge[counter]) + '|' + str(i[1]) + '|' + '{:.5f}'.format(results[i][0]) + '|' + str(len(results[i][1])) +
                          '|' + str(results[i][1]) + '\n' + i[0] + '\n')

            counter += 1

# Save the normalization values in a dataframe
gnorm_df = pd.read_csv(norm_vec_file, index_col=0)

# Change the fasta file into a more workable format
gamma_dic = pu.read_fasta_file(sequences_fasta)

# Calculation of the hydrophobicity + normalization
gamma = g_emul(gamma_dic, gnorm_df)

# Removes peptides depending on the defined cut offs
gamma = cut_offs(gamma)

# Counts each peptides charge
counts = counter(gamma)

# Saves results in viewable file and fasta file for clustering
txt_file(gamma, counts, out_dir)
fasta_file(gamma,counts, out_dir)

