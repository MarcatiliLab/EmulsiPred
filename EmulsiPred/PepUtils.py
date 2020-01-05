import numpy as np
import pandas as pd
import math
import operator


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


# Functions that calculate the hydrophobicity for peptides with random coil
def gamma_emul_calc(seq, Kaa, i):

    # Definition of the hydrophobic and -philic part depending on the cut
    window_l1 = sum([Kaa[seq[n]] for n in range(0 , i)])
    window_l2 = sum([Kaa[seq[n]] for n in range(i ,len(seq))])

    emul = abs(window_l1 - window_l2)

    return(emul)


def g_emul_calcer(seq, Kaa):

    emul_list = []

    # Makes sure the peptide is only cut within its 30% core and has at least 3 amino acids on both sides when cut
    if 3 < ((len(seq) / 2) - 0.15 * len(seq)):
        mini = int((len(seq) / 2)- round(0.15 * len(seq)))
    else:
        mini = 3

    # Calculate the values for each peptide with their selected cut
    for i in range(mini, len(seq)-(mini-1)):
        emul_list.append([gamma_emul_calc(seq, Kaa, i), i])

    return(emul_list)


def g_emul(sequences, norm_df):

    Kaa = hydrophobicity_scale(1)
    g_best = []

    for name in sequences:
        g_scores={}
        seq = sequences[name]

        for i,res in enumerate(seq):

            for j in range(7, 31):
                if i + j <= len(seq):

                    # Calculate the hydrophobicity for the peptide
                    gamma = g_emul_calcer(seq[i:(i+j)], Kaa)
                    gamma = sorted(gamma, key=operator.itemgetter(0), reverse=True)[0:10]

                    # Pick the right values for normalization
                    temp_df = norm_df.loc[(norm_df['Length'] == len(seq[i:(i + j)])) & (norm_df['Cut'] == gamma[0][1])]

                    # Pick cut which gives the highest score (used to represent that peptide) and normalize it
                    g_scores[name, seq[i:(i + j)], gamma[0][1]] = z_normalize(gamma[0][0], temp_df.iloc[0][2], temp_df.iloc[0][3])

        # Make sure it is always only the 10000 best which are saved (makes it run faster)
        g_best = highest_in_list(g_best, g_scores, 10000)

    # Groups repetitive sequences
    switched_gamma = switch_dictionary_gamma(g_best)

    best_emul = pd.DataFrame.from_dict(switched_gamma, orient='index', columns=['score', 'names'])
    best_emul['sequence'] = best_emul.index
    best_emul.reset_index(inplace=True, drop=True)
    best_emul['nr_names'] = best_emul.names.apply(len)
    best_emul['cut'] = best_emul.sequence.apply(lambda x: x[1])
    best_emul['sequence'] = best_emul.sequence.apply(lambda x: x[0])

    return best_emul


# Functions that calculate the hydrophobicity for peptides with a-helices and b-sheets
def accum_emulsifier(emulsify_func, seq, kaa):
    # Function to sum the vectors from alpha and beta
    temp_vec = ([0 ,0])
    window = [emulsify_func(seq[n], n, kaa) for n in range(0, len(seq))]

    for i in window:
        temp_vec[0] += i[0]
        temp_vec[1] += i[1]

    emul_score = math.hypot(temp_vec[0], temp_vec[1])
    return(emul_score)


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


def emul(data_dic, norm_df, type):
    # Main function for calculating emul values for alpha and beta

    threshold = 0.3
    Kaa = hydrophobicity_scale(1)
    best = []

    for pdb_name in data_dic:
        scores, seq, proba = {}, data_dic[pdb_name][0], data_dic[pdb_name][1].split('|')
        proba = [float(i) for i in proba]

        for i,res in enumerate(seq):
            for j in range(7, 30):
                if i + j <= len(seq) and np.mean(proba[i:(i + j)]) >= threshold:
                    # Run equation for the type on the window and sum the vectors
                    emul_val = accum_emulsifier(type, seq[i:(i + j)], Kaa)

                    # Pick the right values for normalization
                    temp_df = norm_df.loc[norm_df['Length'] == len(seq[i:(i + j)])]

                    # Normalize value and save it
                    scores[pdb_name, seq[i:(i + j)]] = z_normalize(emul_val, temp_df.iloc[0][1], temp_df.iloc[0][2])

        # Make sure it is always only the 10000 best which are saved (makes it run faster)
        best = highest_in_list(best, scores, 10000)

    # Groups repetitive sequences
    best_emul = switch_dictionary(best)

    best_emul = pd.DataFrame.from_dict(best_emul, orient='index', columns=['score', 'names'])
    best_emul['sequence'] = best_emul.index
    best_emul.reset_index(inplace=True, drop=True)
    best_emul['nr_names'] = best_emul.names.apply(len)

    return best_emul


# Extra functions
def highest_in_list(knownlist, newlist, number_in_list=100000):

    # The values in newlist is sorted with highest first and the 100 highest is added to knownlist
    knownlist += sorted(newlist.items(), key=operator.itemgetter(1), reverse=True)[0:number_in_list]

    # The values in knownlist is the sorted with highest first and only the 100 highest is saved
    knownlist = sorted(knownlist, key=operator.itemgetter(1), reverse=True)[0:number_in_list]

    return(knownlist)


def switch_dictionary_gamma(thedic):
    # Switches a dictionary around
    lister = {}
    for i in thedic:
        if (i[0][1],i[0][2]) not in lister:
            lister[(i[0][1],i[0][2])] = i[1],[(i[0][0])]
        else:
            for k in lister:
                if (i[0][1],i[0][2]) == k:
                    b = lister[k][1]

                    b = ','.join(b)
                    a = b + ',' + i[0][0]

                    lister[k] = i[1],a.split(','),i[0][2]

    return(lister)


def switch_dictionary(thedic):
# Switches a dictionary around
    lister = {}
    for i in thedic:
        if i[0][1] not in lister:
            lister[i[0][1]] = i[1],[(i[0][0])]
        else:
            for k in lister:
                if i[0][1] == k:
                    b = lister[k][1]

                    b = ','.join(b)
                    a = b + ',' + i[0][0]

                    lister[k] = i[1],a.split(',')

    return(lister)


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