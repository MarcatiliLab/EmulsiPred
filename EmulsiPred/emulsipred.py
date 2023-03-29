import os
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
    
    sequences = read_input(sequences, netsurfp_results=netsurfp_results, peptides=peptides)
    
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


           
def read_input(sequences, netsurfp_results=False, peptides=False):
    
    if not isinstance(sequences, list):
        if os.path.exists(sequences):
            if peptides:
                sequences = pu.read_fasta_file(sequences)
                return list(sequences.values())
            elif not netsurfp_results:
                return pu.read_fasta_file(sequences)
            else:
                return sequences
           
        else:
            sequences = [sequences]
    
    for sequence in sequences:
        for residue in sequence:
            if residue not in 'ABCDEFGHIJKLMNPQRSTVWYXZOU':
                raise ValueError(f'{sequence} is not a valid protein/peptide sequence. Sequences can only contain the following characters: "ABCDEFGHIJKLMNPQRSTVWYXZOU".')

    return sequences        
        