import PepUtils as pu
import argparse
import pandas as pd
import os


class AlphaEmulPred:

    def __init__(self, anorm_file, netsurfp_results, out_dir):
        self.out_dir = out_dir
        # Save the normalization values in a dataframe
        self.norm_df = pd.read_csv(anorm_file, index_col=0)
        # Change the netsurfp results into a more workable format
        self.alpha_dic = pu.get_netsurfp_results(netsurfp_results, 'alpha')
        # Calculation of the hydrophobicity + normalization
        self._predictions = pu.emul(self.alpha_dic, self.norm_df, pu.alpha_emul)
        self._adjusted_predictions = self._predictions

    @property
    def predictions(self):
        # Calculation of the hydrophobicity + normalization
        return self._adjusted_predictions

    def peptide_cutoffs(self, nr_names=4, score=2):
        # Removes peptides depending on the defined cut offs
        self._adjusted_predictions = pu.cut_offs(self._predictions, nr_names, score)

    def save_alpha(self):

        s_df = self._adjusted_predictions
        # Counts each peptides charge
        s_df['charge'] = s_df.sequence.apply(pu.charge_counter)
        # Saves results in a csv format
        s_df.to_csv(os.path.join(self.out_dir, 'a_results.csv'))
        # Saves results in viewable file and fasta file for clustering
        txt_file(self._adjusted_predictions, os.path.join(self.out_dir, 'a_results.txt'))


def txt_file(result_df, out_path):
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


if __name__ == '__main__':

    # python Alpha_emulsifier.py -n ../netsurfp2_results.txt

    # Parse it in
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', dest='netsurfp_results', default='', help='Txt file with netsurfp results')
    parser.add_argument('-v', dest='norm_vals', default=os.path.join('..', 'NormalizationValues', 'a_norm.csv'),
                        help='csv file with mean and standard deviation values')
    parser.add_argument('-o', dest='out_dir', default=os.path.join('..', 'Results'), help='Output directory path')
    parser.add_argument('--n_seq', dest='nr_names', default=1,
                        help='Minimum number of proteins peptides should be a part of.')
    parser.add_argument('--ls', dest='lower_score', default=2., help='Results will only include peptides with a score'
                                                                     'higher than the lowest score.')

    # Define the parsed arguments
    args = parser.parse_args()

    aclass = AlphaEmulPred(args.norm_vals, args.netsurfp_results, args.out_dir)
    aclass.peptide_cutoffs(nr_names=args.nr_names, score=args.lower_score)
    aclass.save_alpha()
