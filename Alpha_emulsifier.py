from __future__ import print_function
import PepUtils as pu
import argparse
import numpy as np
import os
import pandas as pd

# python Alpha_emulsifier.py -n /home/projects/vaccine/people/s132421/Arbejde/Data/netsurfp2_results.txt -v /home/projects/vaccine/people/s132421/Arbejde/Data/NormalizationValues/a_norm.csv -o /home/projects/vaccine/people/s132421/Arbejde/Data/Results

# Parse it in
parser = argparse.ArgumentParser()
parser.add_argument('-n', dest='netsurfp_results', default='', help='Txt file with netsurfp results')
parser.add_argument('-v', dest='norm_vals', default='', help='csv file with mean and standard deviation valuesn')
parser.add_argument('-o', dest='out_dir', default='', help='Output directory path')

# Define the parsed arguments
args = parser.parse_args()
netsurfp_results = args.netsurfp_results
norm_vec_file = args.norm_vals
out_dir = args.out_dir

def cut_offs(results):
    # Removing peptides with less than 2 as a z-score and seen in less than 4 accession numbers
    cut_results = {}
    for i in results:
        if len(results[i][1]) >= 4 and results[i][0] > 2:
            cut_results[i] = results[i]

    return(cut_results)

def counter(alpha):
    # Counts the total charge of each peptide
    total = []
    for i in alpha:
        neg = 0
        pos = 0
        for j in range(len(i)):
            if i[j] == 'D' or i[j] =='E':
                neg += 1
            elif i[j] == 'R' or i[j] == 'K':
                pos += 1

        diff = pos - neg
        total.append(diff)

    return(total)

def txt_file(results, charge, out_dir):
    with open(out_dir+'/a_results.txt', 'w') as outfile:
        outfile.write('-----' + '\t' + 'Charge' + '\t' + 'Kyte-Doolittle' + '\t' + 'Sequence' + '\t\t\t' + 'Diff seqs' + '\n')
        counter = 0
        for i in results:
            outfile.write("ALPHA" + '\t' + str(charge[counter])  +
                          '\t'+ '{:.5f}'.format(results[i][0]) + '\t\t' + "{:<31}".format(i) + '\t' + str(len(results[i][1])) + '\t' + str(results[i][1]) + '\n')

            counter += 1

def fasta_file(results, charge, out_dir):
    with open(out_dir + '/a_results.fsa', 'w') as outfile:
        counter = 0
        for i in results:
            outfile.write(">ALPHA|" + str(charge[counter]) + '|' + '{:.5f}'.format(results[i][0]) + '|' + str(len(results[i][1])) + '|' + str(results[i][1]) + '\n' + i + '\n')

            counter += 1

# Save the normalization values in a dataframe
anorm_df = pd.read_csv(norm_vec_file, index_col=0)

# Change the netsurfp results into a more workable format
alpha_dic = pu.get_netsurfp_results(netsurfp_results,'alpha')

# Calculation of the hydrophobicity + normalization
alpha = pu.emul(alpha_dic,anorm_df,pu.alpha_emul)

# Removes peptides depending on the defined cut offs
alpha = cut_offs(alpha)

# Counts each peptides charge
counts = counter(alpha)

# Saves results in viewable file and fasta file for clustering
txt_file(alpha, counts, out_dir)
fasta_file(alpha, counts, out_dir)


