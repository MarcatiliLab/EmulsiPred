from .predictors import EmulsiPred
import argparse

if __name__ == '__main__':

    # Parse it in
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', dest='sequences', required=True, help='Fasta file with all the sequences.')
    parser.add_argument('-n', dest='netsurfp_results', required=True, help='Txt file with netsurfp results.')
    parser.add_argument('-o', dest='out_dir', default='', help='Directory path for output.')
    parser.add_argument('--nr_seq', dest='nr_seq', default=1,
                        help='Results will only include peptides present in this number of sequences or higher. Default 1.')
    parser.add_argument('--ls', dest='lower_score', default=2.,
                        help='Results will only include peptides with a score higher than this score. Default 2.')

    # Define the parsed arguments
    args = parser.parse_args()

    EmulsiPred(args.sequences, args.netsurfp_results, args.out_dir, args.nr_seq, args.lower_score)