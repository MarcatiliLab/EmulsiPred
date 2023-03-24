import numpy as np
import pandas as pd
import math
import operator
from glob import glob
import os

import pkg_resources


def hydrophobicity_scale(hydro):

    # Hydrophobicity scale for Kyte-Doolittle:
    if hydro == 1:
        Kaa = {'A': 1.80, 'B': 0, 'C': 2.50, 'D': -3.50, 'E': -3.50, 'F': 2.80, 'G': -0.40, 'H': -3.20, 'I': 4.50,
               'J': 0, 'K': -3.90, 'L': 3.80, 'M': 1.90, 'N': -3.50, 'P': -1.60, 'Q': -3.50, 'R': -4.50, 'S': -0.80,
               'T': -0.70, 'V': 4.20, 'W': -0.90, 'Y': -1.30, 'X': 0, 'Z': 0, 'O': 0, 'U': 0}
    # Hydrophobicity scale for Hopp-Woods:
    if hydro == 2:
        Kaa = {'A': -0.50, 'B': 0, 'C': -1.00, 'D': 3.00, 'E': 3.00, 'F': -2.50, 'G': 0.00, 'H': -0.50, 'I': -1.80,
               'J': 0, 'K': 3.00, 'L': -1.80, 'M': -1.30, 'N': 0.20, 'P': 0.00, 'Q': 0.20, 'R': 3.00, 'S': 0.30,
               'T': -0.40, 'V': -1.50, 'W': -3.40, 'Y': -2.30, 'X': 0, 'Z': 0, 'O': 0, 'U': 0}

    return Kaa


def z_normalize(score, mean, std):
    # Normalizes with the help of z-score
    normalized = (score-mean)/std
    return(normalized)


# Read inputs (sequences in fasta file and netsurfp results)
def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))


def read_fasta_file(fasta_file):
    gamma_dic = {}
    with open(fasta_file) as sequences:
        for name, seq in read_fasta(sequences):
            gamma_dic[name[1:]] = seq

    return(gamma_dic)


def get_netsurfp_results(netsurfp_results, type):

    resi_number = 0
    seqs = ''
    score = []
    dic_of_vals = {}

    for line in open(netsurfp_results, 'r'):
        line = line.split(' ')
        line = list(filter(None, line))

        if line[0] == 'E' or line[0] == 'B':
            # make new list for next sequence
            if (int(line[3]) - resi_number) < 0:
                seqs = ''
                score = []
                resi_number = int(line[3])


            seqs += line[1]
            # Save the score depending on which conformation you are looking at
            if type == 'alpha':
                score.append(line[7])
            elif type == 'beta':
                score.append(line[8])

            dic_of_vals[line[2]] = seqs, '|'.join(score)
            resi_number += 1

    return(dic_of_vals)


def get_netsurfp_csv(adir, sstype='alpha'):
    
    dic_of_vals = {}

    for file in glob(os.path.join(adir, "*.csv")):
        outcsv = pd.read_csv(file, index_col=0)
        outcsv = outcsv.astype(str)
        
        if sstype=='alpha':
            sstype= "q3(h)"
        elif sstype=='beta':
            sstype="q3(e)"
        
        score = '|'.join(list(outcsv[sstype].values))
        seqs = ''.join(list(outcsv["residue"].values))

        dic_of_vals[os.path.basename(file[:-4])] = seqs, score
    
    return dic_of_vals


def get_peptides(sequences, gamma_emul = False, threshold = 0.3):
    
    if gamma_emul:
        for name, seq in sequences.items():
            for i,res in enumerate(seq):
                for j in range(7, 31):
                    if i + j <= len(seq): yield name, seq[i:(i+j)]
    
    else:
        for name, (seq, secondary_structure) in sequences.items():
            
            secondary_structure = [float(i) for i in secondary_structure.split('|')]

            for i,res in enumerate(seq):
                for j in range(7, 31):
                    if i + j <= len(seq) and np.mean(secondary_structure[i:(i + j)]) >= threshold: 
                        yield name, seq[i:(i+j)]
                    

# Functions that calculate the hydrophobicity for peptides with random coil
def calculate_gamma_emul_single_cut(seq, Kaa, i):

    # Definition of the hydrophobic and -philic part depending on the cut
    window_l1 = sum([Kaa[seq[n]] for n in range(0 , i)])
    window_l2 = sum([Kaa[seq[n]] for n in range(i ,len(seq))])

    emul = abs(window_l1 - window_l2)

    return emul


def calculate_gamma_emul(seq, Kaa, gamma_norm):

    final_emul = [-100, 0]

    # Makes sure the peptide is only cut within its 30% core and has at least 3 amino acids on both sides when cut
    if 3 < ((len(seq) / 2) - 0.15 * len(seq)):
        mini = int((len(seq) / 2)- round(0.15 * len(seq)))
    else:
        mini = 3

    # Calculate the values for each peptide with their selected cut
    for i in range(mini, len(seq)-(mini-1)):
        gamma_emul = calculate_gamma_emul_single_cut(seq, Kaa, i)
        
        mean, std = gamma_norm[len(seq)][i]
        gamma_emul = z_normalize(gamma_emul, mean, std)

        if final_emul[0] < gamma_emul: final_emul = [gamma_emul, i]
        
    return final_emul

def gamma_norm_df_to_norm_dict(norm_df):
    norm_dict = {}
    for length, grp in norm_df.groupby(by='Length'):
        tmp_dict = {}
        for _, cut, mean, std in grp.values:
            tmp_dict[int(cut)] = [mean,std]

        norm_dict[length] = tmp_dict

    return norm_dict

def many_calculate_gamma_emul(sequences, norm_df):

    Kaa = hydrophobicity_scale(1)
    gamma_norm = gamma_norm_df_to_norm_dict(norm_df)

    results = []
    for name, seq in get_peptides(sequences, gamma_emul = True):
        results.append([name, seq, *calculate_gamma_emul(seq, Kaa, gamma_norm)])

    results = pd.DataFrame(results, columns = ['name', 'sequence', 'score', 'cut']).sort_values(by='score', ascending=False).iloc[:10000] 

    groupped_results = []
    for seq, grp in results.groupby(by='sequence'):
        #, seq, grp.score[0], grp.cut[0])
        groupped_results.append([grp.score.values[0], grp.name.values, seq, len(grp.name.values), grp.cut.values[0]])

    groupped_results = pd.DataFrame(groupped_results, columns = ['score', 'names', 'sequence', 'nr_names', 'cut']) 

    return groupped_results.sort_values(by='score', ascending=False).reset_index(drop=True)


# Functions that calculate the hydrophobicity for peptides with a-helices and b-sheets
def calculate_emul(seq, kaa, emulsify_func, emul_norm):

    # Function to sum the vectors from alpha and beta
    temp_vec = ([0 ,0])
    window = [emulsify_func(seq[n], n, kaa) for n in range(0, len(seq))]

    for i in window:
        temp_vec[0] += i[0]
        temp_vec[1] += i[1]

    emul_score = math.hypot(temp_vec[0], temp_vec[1])
    
    mean, std = emul_norm[len(seq)]
    return z_normalize(emul_score, mean, std)


def alpha_emul(amino, n, Kaa):
    # Find the hydrophobicity value for each amino acid
    n += 1
    K = Kaa[amino]
    vector = K * np.array([math.cos(n * (( 5 /9 ) *math.pi)), math.sin(n * (( 5 /9 ) *math.pi))])
    return vector


def beta_emul(amino, n, Kaa):
    # Find the hydrophobicity value for each amino acid
    n += 1
    K = Kaa[amino]
    vector = K * np.array([math.cos(n * (math.pi)), math.sin(n * (math.pi))])
    return vector


def many_calculate_emul(sequences, norm_df, emul_fn):

    Kaa = hydrophobicity_scale(1)
    emul_norm = {}
    for length, mean, std in norm_df.values: emul_norm[length] = [mean, std]

    results = []
    for name, seq in get_peptides(sequences):
        results.append([name, seq, calculate_emul(seq, Kaa, emul_fn, emul_norm)])

    results = pd.DataFrame(results, columns = ['name', 'sequence', 'score']).sort_values(by='score', ascending=False).iloc[:10000] 

    groupped_results = []
    for seq, grp in results.groupby(by='sequence'):
        groupped_results.append([grp.score.values[0], grp.name.values, seq, len(grp.name.values)])

    groupped_results = pd.DataFrame(groupped_results, columns = ['score', 'names', 'sequence', 'nr_names']) 

    return groupped_results.sort_values(by='score', ascending=False).reset_index(drop=True)


def peptide_predicter(peptides):
    
    Kaa = hydrophobicity_scale(1)
    
    anormdf = pd.read_csv(pkg_resources.resource_filename(__name__, os.path.join('NormalizationValues', 'a_norm.csv')), index_col=0)
    bnormdf = pd.read_csv(pkg_resources.resource_filename(__name__, os.path.join('NormalizationValues', 'b_norm.csv')), index_col=0)
    gnormdf = pd.read_csv(pkg_resources.resource_filename(__name__, os.path.join('NormalizationValues', 'g_norm.csv')), index_col=0)
    
    alpha_norm = {}
    for length, mean, std in anormdf.values: alpha_norm[length] = [mean, std]
    
    beta_norm = {}
    for length, mean, std in bnormdf.values: beta_norm[length] = [mean, std]
    
    gamma_norm = gamma_norm_df_to_norm_dict(gnormdf)
    
    
    peptides['alpha'] = peptides.seq.apply(lambda x: calculate_emul(x, Kaa, alpha_emul, alpha_norm))
    peptides['beta'] = peptides.seq.apply(lambda x: calculate_emul(x, Kaa, beta_emul, beta_norm))
    peptides[['gamma', 'cut']] = peptides.seq.apply(lambda x: pd.Series(calculate_gamma_emul(x, Kaa, gamma_norm)))
    
    return peptides
  

def cut_offs(results, nr_seq=4, score=2):
    # Removing peptides with less than 2 as a z-score and seen in less than 4 accession numbers

    bdf = results.copy()

    bdf = bdf[bdf['nr_names'] >= nr_seq]
    bdf = bdf[bdf['score'] >= score]

    return bdf


def charge_counter(seq):
    neg = 0
    pos = 0
    for res in range(len(seq)):
        if seq[res] == 'D' or seq[res] =='E':
            neg += 1
        elif seq[res] == 'R' or seq[res] == 'K':
            pos += 1

    return pos - neg


# Save functions
def a_txt_file(result_df, out_path):
    with open(out_path, 'w') as outfile:
        outfile.write('-----' + '\t' +
                      'Charge' + '\t' +
                      'Kyte-Doolittle' + '\t' +
                      'Sequence' + '\t\t\t' +
                      'Diff seqs' + '\n')
        for row in result_df.itertuples(index=True, name='Pandas'):
            outfile.write("ALPHA" + '\t' +
                          "{:<5}".format(str(getattr(row, "charge"))) + '\t' +
                          '{:.5f}'.format(getattr(row, 'score')) + '\t\t\t' +
                          "{:<31}".format(getattr(row, 'sequence')) + '\t' +
                          str(getattr(row, 'nr_names')) + '\t' +
                          str(getattr(row, 'names')) + '\n')


def b_txt_file(result_df, out_path):
    with open(out_path, 'w') as outfile:
        outfile.write('-----' + '\t' +
                      'Charge' + '\t' +
                      'Kyte-Doolittle' + '\t' +
                      'Sequence' + '\t\t\t' +
                      'Diff seqs' + '\n')
        for row in result_df.itertuples(index=True, name='Pandas'):
            outfile.write("BETA" + '\t' +
                          "{:<5}".format(str(getattr(row, "charge"))) + '\t' +
                          '{:.5f}'.format(getattr(row, 'score')) + '\t\t\t' +
                          "{:<31}".format(getattr(row, 'sequence')) + '\t' +
                          str(getattr(row, 'nr_names')) + '\t' +
                          str(getattr(row, 'names')) + '\n')


def g_txt_file(result_df, out_path):
    # Writes the results in a nice viewable file
    with open(out_path, 'w') as outfile:
        outfile.write('-----' + '\t' +
                      'Charge' + '\t' +
                      'Cut' '\t\t' +
                      'Kyte-Doolittle' + '\t' +
                      'Sequence' + '\t\t\t' +
                      'Diff seqs' + '\n')

        for row in result_df.itertuples(index=True, name='Pandas'):
            outfile.write("GAMMA" + '\t' +
                          "{:<5}".format(str(getattr(row, "charge"))) + '\t' +
                          '{}'.format(getattr(row, 'cut')) + '\t\t' +
                          '{:.5f}'.format(getattr(row, 'score')) + '\t\t\t' +
                          "{:<31}".format(getattr(row, 'sequence')) + '\t' +
                          str(getattr(row, 'nr_names')) + '\t' +
                          str(getattr(row, 'names')) + '\n')