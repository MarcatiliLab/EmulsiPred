import PepUtils as pu
import argparse
import pandas as pd
import os

class GammaEmulPred:

    def __init__(self, gnorm_file, sequence_fsa, out_dir):
        self.out_dir = out_dir
        # Save the normalization values in a dataframe
        self.norm_df = pd.read_csv(gnorm_file, index_col=0)
        # Change the netsurfp results into a more workable format
        self.gamma_dic = pu.read_fasta_file(sequence_fsa)
        # Calculation of the hydrophobicity + normalization
        self._predictions = pu.g_emul(self.gamma_dic, self.norm_df)
        self._adjusted_predictions = self._predictions

    @property
    def predictions(self):
        # Calculation of the hydrophobicity + normalization
        return self._adjusted_predictions

    def peptide_cutoffs(self, nr_names=4, score=2.):
        # Removes peptides depending on the defined cut offs
        self._adjusted_predictions = pu.cut_offs(self._predictions, nr_names, score)

    def save_gamma(self):
        s_df = self._adjusted_predictions
        # Counts each peptides charge
        s_df['charge'] = s_df.sequence.apply(pu.charge_counter)
        # Saves results in a csv format
        s_df.to_csv(os.path.join(self.out_dir, 'g_results.csv'))
        # Saves results in viewable file and fasta file for clustering
        txt_file(self._adjusted_predictions, os.path.join(self.out_dir, 'g_results.txt'))


def txt_file(result_df, out_path):
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



if __name__ == '__main__':

    # python Gamma_emulsifier.py -s ../PeptideSequences.fsa

    # Parse it in
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', dest='sequences', default='', help='Fasta file with all the sequences')
    parser.add_argument('-v', dest='norm_vals', default=os.path.join('..', 'NormalizationValues', 'g_norm.csv'),
                        help='csv file with mean and standard deviation values')
    parser.add_argument('-o', dest='out_dir', default=os.path.join('..', 'Results'), help='Output directory path')
    parser.add_argument('--n_seq', dest='nr_names', default=1,
                        help='Minimum number of proteins peptides should be a part of.')
    parser.add_argument('--ls', dest='lower_score', default=2., help='Results will only include peptides with a score'
                                                                  'higher than the lowest score.')

    # Define the parsed arguments
    args = parser.parse_args()

    gclass = GammaEmulPred(args.norm_vals, args.sequences, args.out_dir)
    gclass.peptide_cutoffs(nr_names=int(args.nr_names), score=float(args.lower_score))
    gclass.save_gamma()

