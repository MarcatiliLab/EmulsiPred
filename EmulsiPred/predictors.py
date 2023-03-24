import EmulsiPred.utils as pu
import pandas as pd
import os
import pkg_resources


class AlphaEmulPred:

    def __init__(self, sequences, out_dir, netsurfp_results):
        self.out_dir = out_dir
        # Save the normalization values in a dataframe
        self.norm_df = pd.read_csv(pkg_resources.resource_filename(
            __name__, os.path.join('NormalizationValues', 'a_norm.csv')), index_col=0)
        # Change the netsurfp results into a more workable format
        
        if netsurfp_results:
            if os.path.isdir(sequences):
                self.alpha_dic = pu.get_netsurfp_csv(sequences, 'alpha')
            else:
                self.alpha_dic = pu.get_netsurfp_results(sequences, 'alpha')
        else:
            self.alpha_dic = sequences.copy()
            for key, value in self.alpha_dic.items():
                self.alpha_dic[key] = value, '|'.join(['1.000' for _ in value])
        
        # Calculation of the hydrophobicity + normalization
        self._predictions = pu.many_calculate_emul(self.alpha_dic, self.norm_df, pu.alpha_emul)
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
        s_df.to_csv(os.path.join(self.out_dir, 'a_results.csv'), index=0)
        # Saves results in viewable file and fasta file for clustering
        pu.a_txt_file(self._adjusted_predictions, os.path.join(self.out_dir, 'a_results.txt'))


class BetaEmulPred:

    def __init__(self, sequences, out_dir, netsurfp_results):
        self.out_dir = out_dir
        # Save the normalization values in a dataframe
        self.norm_df = pd.read_csv(pkg_resources.resource_filename(
            __name__, os.path.join('NormalizationValues', 'b_norm.csv')), index_col=0)
        # Change the netsurfp results into a more workable format           
        if netsurfp_results:
            if os.path.isdir(sequences):
                self.alpha_dic = pu.get_netsurfp_csv(sequences, 'beta')
            else:
                self.alpha_dic = pu.get_netsurfp_results(sequences, 'beta')
        else:
            self.alpha_dic = sequences.copy()
            for key, value in self.alpha_dic.items():
                self.alpha_dic[key] = value, '|'.join(['1.000' for _ in value])
            
        # Calculation of the hydrophobicity + normalization
        self._predictions = pu.many_calculate_emul(self.alpha_dic, self.norm_df, pu.beta_emul)
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
        s_df.to_csv(os.path.join(self.out_dir, 'b_results.csv'), index=0)
        # Saves results in viewable file and fasta file for clustering
        pu.b_txt_file(self._adjusted_predictions, os.path.join(self.out_dir, 'b_results.txt'))


class GammaEmulPred:

    def __init__(self, sequences, out_dir, netsurfp_results):
        self.out_dir = out_dir
        # Save the normalization values in a dataframe
        self.norm_df = pd.read_csv(pkg_resources.resource_filename(
            __name__, os.path.join('NormalizationValues', 'g_norm.csv')), index_col=0)
        # Change the netsurfp results into a more workable format
        if netsurfp_results:
            if os.path.isdir(sequences):
                self.gamma_dic = pu.get_netsurfp_csv(sequences, 'beta')
                
            else:
                self.gamma_dic = pu.get_netsurfp_results(sequences, 'beta')

            for key, value in self.gamma_dic.items():
                self.gamma_dic[key] = value[0]
                
        else:
            self.gamma_dic = sequences.copy()
        
        # Calculation of the hydrophobicity + normalization
        self._predictions = pu.many_calculate_gamma_emul(self.gamma_dic, self.norm_df)
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
        s_df.to_csv(os.path.join(self.out_dir, 'g_results.csv'), index=0)
        # Saves results in viewable file and fasta file for clustering
        pu.g_txt_file(self._adjusted_predictions, os.path.join(self.out_dir, 'g_results.txt'))