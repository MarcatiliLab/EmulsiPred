import argparse
from distutils.util import strtobool

from EmulsiPred import EmulsiPred

if __name__ == '__main__':

    # Parse it in
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', dest='sequences',  type=str, required=True, help='Fasta file with all the sequences or txt file with netsurfp results.')
    parser.add_argument('-n', dest='netsurfp_results',  type=strtobool, default=False, help='Whether input is netsurfp results')
    parser.add_argument('-p', dest='peptides',  type=strtobool, default=False, help='Whether input is peptides (this will omit splitting input into smaller peptides)')
    parser.add_argument('-o', dest='out_folder',  type=str, default='', help='Directory path for output.')
    parser.add_argument('--nr_seq', dest='nr_seq',  type=int, default=1,
                        help='Results will only include peptides present in this number of sequences or higher. Default 1.')
    parser.add_argument('--ls', dest='lower_score',  type=int, default=2.,
                        help='Results will only include peptides with a score higher than this score. Default 2.')

    # Define the parsed arguments
    args = parser.parse_args()

    EmulsiPred(args.sequences, args.netsurfp_results, args.peptides, args.out_folder, args.nr_seq, args.lower_score)