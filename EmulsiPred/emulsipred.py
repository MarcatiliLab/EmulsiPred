import argparse
import os
from distutils.util import strtobool

import pandas as pd

from EmulsiPred import AlphaEmulPred, BetaEmulPred, GammaEmulPred
import EmulsiPred.utils as pu


def EmulsiPred(sequences, netsurfp_results=False, peptides=False, out_folder='', seen_in_N_seqs=1, lowest_score=2):
    """
    sequences: Either a file (fasta or netsurfp), list of peptides/sequences or a string.
    netsurfp_results: True if file is netsurfp, otherwise False.
    peptides: True if input is to be treated as peptides, otherwise False. If treated as a peptide, predictions will only be made for that specific peptide and not windows of the peptide as well (as done for sequences). A peptide is defined as 7-30 aa's in length (peptides outside this length will be removed).
    out_folder: Specific folder to save data in.
    seen_in_N_seqs: Sequence only argument. Keep only results seen in at least N number of sequences.
    lowest_score: Sequence only argument. Remove results with scores lower than this value.
    """
    os.makedirs(out_folder, exist_ok=True)
    
    if not isinstance(sequences, list):
        if os.path.isfile(sequences) and (netsurfp_results==False) and (peptides==False):
            sequences = pu.read_fasta_file(sequences)
            
        elif os.path.isfile(sequences) and (peptides==True):
            sequences = pu.read_fasta_file(sequences)
            sequences = list(sequences.values())
        
        elif (os.path.isfile(sequences) or os.path.isdir(sequences)) and (netsurfp_results==True):
            pass
        
        else:
            sequences = [sequences]    
    
    if peptides==True:
        sequences = pd.DataFrame(sequences, columns=['seq']).query("seq.str.len()<30", engine='python')
        results = pu.peptide_predicter(sequences)
        
        results['charge'] = results.seq.apply(pu.charge_counter)
        results.to_csv(os.path.join(out_folder, "emul_results.csv"), index=False)
    
    else:
        a_class = AlphaEmulPred(sequences, out_folder, netsurfp_results)
        a_class.peptide_cutoffs(nr_seq=int(seen_in_N_seqs), score=float(lowest_score))
        a_class.save_alpha()

        b_class = BetaEmulPred(sequences, out_folder, netsurfp_results)
        b_class.peptide_cutoffs(nr_seq=int(seen_in_N_seqs), score=float(lowest_score))
        b_class.save_beta()

        g_class = GammaEmulPred(sequences, out_folder, netsurfp_results)
        g_class.peptide_cutoffs(nr_seq=int(seen_in_N_seqs), score=float(lowest_score))
        g_class.save_gamma()


if __name__ == '__main__':

    # Parse it in
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', dest='sequences',  type=str, required=True, help='Fasta file with all the sequences or txt file with netsurfp results.')
    parser.add_argument('-n', dest='netsurfp_results',  type=strtobool, default=False, help='Whether input is netsurfp results')
    parser.add_argument('-p', dest='peptides',  type=strtobool, default=False, help='Whether input is peptides (this will omit splitting input into smaller peptides)')
    parser.add_argument('-o', dest='out_dir',  type=str, default='', help='Directory path for output.')
    parser.add_argument('--nr_seq', dest='nr_seq',  type=int, default=1,
                        help='Results will only include peptides present in this number of sequences or higher. Default 1.')
    parser.add_argument('--ls', dest='lower_score',  type=int, default=2.,
                        help='Results will only include peptides with a score higher than this score. Default 2.')

    # Define the parsed arguments
    args = parser.parse_args()

    EmulsiPred(args.sequences, args.netsurfp_results, args.peptides, args.out_dir, args.nr_seq, args.lower_score)