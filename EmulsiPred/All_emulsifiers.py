from .Alpha_emulsifier import AlphaEmulPred
from .Beta_emulsifier import BetaEmulPred
from .Gamma_emulsifier import GammaEmulPred
import argparse
import os

def EmulsiPred(sequences, netsurfp_results, norm_vals, out_dir, nr_names, lower_score):

    aclass = AlphaEmulPred(os.path.join(norm_vals, 'a_norm.csv'), netsurfp_results, out_dir)
    aclass.peptide_cutoffs(nr_names=int(nr_names), score=float(lower_score))
    aclass.save_alpha()

    bclass = BetaEmulPred(os.path.join(norm_vals, 'b_norm.csv'), netsurfp_results, out_dir)
    bclass.peptide_cutoffs(nr_names=int(nr_names), score=float(lower_score))
    bclass.save_beta()

    gclass = GammaEmulPred(os.path.join(norm_vals, 'g_norm.csv'), sequences, out_dir)
    gclass.peptide_cutoffs(nr_names=int(nr_names), score=float(lower_score))
    gclass.save_gamma()


if __name__ == '__main__':

    # python All_emulsifiers.py -n ../netsurfp2_results.txt -v ../a_norm.csv -o ../Results

    # Parse it in
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', dest='sequences', default='', help='Fasta file with all the sequences')
    parser.add_argument('-n', dest='netsurfp_results', default='', help='Txt file with netsurfp results')
    parser.add_argument('-v', dest='norm_vals', default=os.path.join('..', 'NormalizationValues'),
                        help='Folder with csv files with mean and standard deviation values')
    parser.add_argument('-o', dest='out_dir', default=os.path.join('..', 'Results'), help='Output directory path')
    parser.add_argument('--n_seq', dest='nr_names', default=1,
                        help='Minimum number of proteins peptides should be a part of.')
    parser.add_argument('--ls', dest='lower_score', default=2., help='Results will only include peptides with a score'
                                                                     'higher than the lowest score.')

    # Define the parsed arguments
    args = parser.parse_args()

    EmulsiPred(args.sequences, args.netsurfp_results, args.norm_vals, args.out_dir, args.nr_names, args.lower_score)