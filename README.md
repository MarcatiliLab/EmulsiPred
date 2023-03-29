# EmulsiPred
Tool for prediction of emulsifying peptides. EmulsiPred predicts the emulsifying property of either a single peptide or for any peptide within a protein sequences. Three emulsifying scores are calculated for each peptide as described by [García-Moreno P.J. et al., 2020](doi.org/10.1038/s41598-019-57229-6), with a peptide defined as a sequence of 7-30 amino acids.

EmulsiPred takes as input a fasta file or a NetSurfP (2 or 3) result file. The NetSurfP-2 file should be in the NetSurfP-1 Format (retrieved when clicking 'Export All' in the upper right side of NetSurfP's 'Server Output' window). For a fasta file with protein sequences, EmulsiPred will return scores for each peptide found within the protein sequences. If given a NetSurfP result file, EmulsiPred will only return the alpha and beta scores for peptides present in either an alpha helix or beta sheet, predicted by NetSurfP. 


#### Prerequisites and installation

The package can either be installed with pip or from github. 
In both cases, python-3.9 or higher needs to be installed in your
environment. Additionally, it is 
recommended to install the package in a new environment.

The following commands are run in the command line.

1: Set up a new environment.
~~~.sh  
    python3 -m venv EmulsiPred_env
~~~
2: Enter (activate) the environment.
~~~.sh
    source EmulsiPred_env/bin/activate
~~~
3a: Install EmulsiPred within the activated environment with pip.
~~~.sh
    pip install EmulsiPred
~~~
    
3b: Install EmulsiPred by installing from github with pip.

~~~.sh
    pip install "git+https://github.com/MarcatiliLab/EmulsiPred.git"
~~~ 

After either running 3a or 3b, EmulsiPred is installed within the
activated environment (in our case EmulsiPred_env).

---
#### Running EmulsiPred

After installation, EmulsiPred can be run from the terminal or
within a python script.

As mentioned above, EmulsiPred requires a fasta file containing the protein sequences to check for emulsifiers or a NetSurfP file containing secondary structure information of each sequences.

Additionally, there are also five additional parameters. 
1) -n (netsurfp_results): Whether the input is a NetSurfP file (default is False)
2) -p (peptides): Whether the input are peptides and therefore shouldn't be cleaved into peptides (default is False)
3) -o (out_dir): Output directory (default is the current directory).
4) --nr_seq (nr_seq): Results will only include peptides present in this number of sequences or higher (default 1).
5) --ls (lower_score): Results will only include peptides with a score higher than this score (default 2).  

EmulsiPred can be run directly in the terminal with the following
command.
~~~.sh
    python -m EmulsiPred -s path/to/sequence.fsa -n False -p False -o path/to/out_dir --nr_seq 1 --ls 2
~~~ 
Furthermore, it can be imported and run in a python script.

~~~~~~~~~~~~~~~~~~~~~python
import EmulsiPred as ep

ep.EmulsiPred(sequences='path/to/sequence.fsa', netsurfp_results=False, peptides=False, out_dir='path/to/out_dir', nr_seq=1, lower_score=2)
~~~~~~~~~~~~~~~~~~~~~

#### Interpretation of predictions

The predicted values are a relative ordering of the peptides by chance of being an emulsifier. 
In other words, a higher score implies a higher chance of being an emulsifier. 