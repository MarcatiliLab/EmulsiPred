import EmulsiPred.PepUtils as pu
import pandas as pd
import os
import pkg_resources


def EmulsiPred(sequences, netsurfp_results, out_dir='', nr_seq='1', lower_score='2.'):

    a_class = AlphaEmulPred(netsurfp_results, out_dir)
    a_class.peptide_cutoffs(nr_seq=int(nr_seq), score=float(lower_score))
    a_class.save_alpha()

    b_class = BetaEmulPred(netsurfp_results, out_dir)
    b_class.peptide_cutoffs(nr_seq=int(nr_seq), score=float(lower_score))
    b_class.save_beta()

    g_class = GammaEmulPred(sequences, out_dir)
    g_class.peptide_cutoffs(nr_seq=int(nr_seq), score=float(lower_score))
    g_class.save_gamma()


class AlphaEmulPred:

    def __init__(self, netsurfp_results, out_dir):
        self.out_dir = out_dir
        # Save the normalization values in a dataframe
        self.norm_df = pd.read_csv(pkg_resources.resource_filename(
            __name__, os.path.join('NormalizationValues', 'a_norm.csv')), index_col=0)
        # Change the netsurfp results into a more workable format
        self.alpha_dic = pu.get_netsurfp_results(netsurfp_results, 'alpha')
        # Calculation of the hydrophobicity + normalization
        self._predictions = pu.emul(self.alpha_dic, self.norm_df, pu.alpha_emul)
        self._adjusted_predictions = self._predictions

    @property
    def predictions(self):
        # Calculation of the hydrophobicity + normalization
        return self._adjusted_predictions

    def peptide_cutoffs(self, nr_seq=4, score=2.):
        # Removes peptides depending on the defined cut offs
        self._adjusted_predictions = pu.cut_offs(self._predictions, nr_seq, score)

    def save_alpha(self):

        s_df = self._adjusted_predictions
        # Counts each peptides charge
        s_df['charge'] = s_df.sequence.apply(pu.charge_counter)
        # Saves results in a csv format
        s_df.to_csv(os.path.join(self.out_dir, 'a_results.csv'))
        # Saves results in viewable file and fasta file for clustering
        pu.a_txt_file(self._adjusted_predictions, os.path.join(self.out_dir, 'a_results.txt'))


class BetaEmulPred:

    def __init__(self, netsurfp_results, out_dir):
        self.out_dir = out_dir
        # Save the normalization values in a dataframe
        self.norm_df = pd.read_csv(pkg_resources.resource_filename(
            __name__, os.path.join('NormalizationValues', 'b_norm.csv')), index_col=0)
        # Change the netsurfp results into a more workable format
        self.alpha_dic = pu.get_netsurfp_results(netsurfp_results, 'beta')
        # Calculation of the hydrophobicity + normalization
        self._predictions = pu.emul(self.alpha_dic, self.norm_df, pu.alpha_emul)
        self._adjusted_predictions = self._predictions

    @property
    def predictions(self):
        # Calculation of the hydrophobicity + normalization
        return self._adjusted_predictions

    def peptide_cutoffs(self, nr_seq=4, score=2.):
        # Removes peptides depending on the defined cut offs
        self._adjusted_predictions = pu.cut_offs(self._predictions, nr_seq, score)

    def save_beta(self):
        s_df = self._adjusted_predictions
        # Counts each peptides charge
        s_df['charge'] = s_df.sequence.apply(pu.charge_counter)
        # Saves results in a csv format
        s_df.to_csv(os.path.join(self.out_dir, 'b_results.csv'))
        # Saves results in viewable file and fasta file for clustering
        pu.b_txt_file(self._adjusted_predictions, os.path.join(self.out_dir, 'b_results.txt'))


class GammaEmulPred:

    def __init__(self, sequence_fsa, out_dir):
        self.out_dir = out_dir
        # Save the normalization values in a dataframe
        self.norm_df = pd.read_csv(pkg_resources.resource_filename(
            __name__, os.path.join('NormalizationValues', 'g_norm.csv')), index_col=0)
        # Change the netsurfp results into a more workable format
        self.gamma_dic = pu.read_fasta_file(sequence_fsa)
        # Calculation of the hydrophobicity + normalization
        self._predictions = pu.g_emul(self.gamma_dic, self.norm_df)
        self._adjusted_predictions = self._predictions

    @property
    def predictions(self):
        # Calculation of the hydrophobicity + normalization
        return self._adjusted_predictions

    def peptide_cutoffs(self, nr_seq=4, score=2.):
        # Removes peptides depending on the defined cut offs
        self._adjusted_predictions = pu.cut_offs(self._predictions, nr_seq, score)

    def save_gamma(self):
        s_df = self._adjusted_predictions
        # Counts each peptides charge
        s_df['charge'] = s_df.sequence.apply(pu.charge_counter)
        # Saves results in a csv format
        s_df.to_csv(os.path.join(self.out_dir, 'g_results.csv'))
        # Saves results in viewable file and fasta file for clustering
        pu.g_txt_file(self._adjusted_predictions, os.path.join(self.out_dir, 'g_results.txt'))