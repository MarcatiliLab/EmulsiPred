import EmulsiPred.PepUtils as pu
import pandas as pd
import os
import pkg_resources


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
        anormdf = pd.read_csv(pkg_resources.resource_filename(__name__, os.path.join('NormalizationValues', 'a_norm.csv')), index_col=0)
        bnormdf = pd.read_csv(pkg_resources.resource_filename(__name__, os.path.join('NormalizationValues', 'b_norm.csv')), index_col=0)
        gnormdf = pd.read_csv(pkg_resources.resource_filename(__name__, os.path.join('NormalizationValues', 'g_norm.csv')), index_col=0)
        results = pu.peptide_predicter(sequences, anormdf, bnormdf, gnormdf)
        
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
        self._predictions = pu.emul(self.alpha_dic, self.norm_df, pu.beta_emul)
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
        s_df.to_csv(os.path.join(self.out_dir, 'g_results.csv'), index=0)
        # Saves results in viewable file and fasta file for clustering
        pu.g_txt_file(self._adjusted_predictions, os.path.join(self.out_dir, 'g_results.txt'))