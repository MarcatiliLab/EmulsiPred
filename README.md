# EmulsiPred
Prediction of Emulsifying Peptides, based on protein sequences (in fasta format) and
their corresponding results from NetSurfP-2 (http://www.cbs.dtu.dk/services/NetSurfP/).
The NetSurfP-2 file should be in the NetSurfP-1 Format (retrieved when clicking 'Export All'
in the upper right side of NetSurfP's 'Server Output' window).


#### Prerequisites and installation

The package can either be cloned from github and installed 
locally or installed with pip. In both cases, python-3.6 or 
higher needs to be installed on your PC. Additionally, it is 
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

After either running 3a or 3b EmulsiPred is installed within the
activated environment (in our case EmulsiPred_env).

---
#### Running EmulsiPred

After installation, EmulsiPred can be run from the terminal or
within a python script.

As mentioned above, EmulsiPred requires 2 inputs.
1) A fasta file containing the protein sequences to check for emulsifiers (termed sequence.fsa) or a NetSurfP file containing secondary structure information of the sequences in sequence.fsa (termed netsurfp.txt)  
2) Whether the input is a netsurfp file
3) Whether the input are peptides (and therefore shouldn't be cleaved into peptides)

Additionally, there are also 3 variable parameters. 
1) o (out_dir): Output directory (default is the current directory).
2) nr_seq: Results will only include peptides present in this number of sequences or higher (default 1).
3) ls (lower score): Results will only include peptides with a score higher than this score (default 2).  

EmulsiPred can be run directly in the terminal with the following
command.
~~~.sh
    python -m EmulsiPred -s path/to/sequence.fsa -n False -p False -o path/to/out_dir --nr_seq 1 --ls 2
~~~ 
Furthermore, it can be imported and run in a python script.

~~~~~~~~~~~~~~~~~~~~~python
import EmulsiPred as ep

ep.EmulsiPred(sequences='path/to/sequence.fsa', netsurfp_results=False, peptide=False, out_dir='path/to/out_dir', nr_seq=1, lower_score=2)
~~~~~~~~~~~~~~~~~~~~~

#### Interpretation of predictions

The predicted values are a relative ordering 
of the peptides by chance of being an emulsifier. 
In other words, a higher score implies a higher chance 
of being an emulsifier. 