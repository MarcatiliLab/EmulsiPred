# EmulsiPred
Prediction of Emulsifying Peptides


#### Prerequisites and installation

The package can either be cloned from github and installed 
locally or installed with pip. In both cases, python-3.6 or 
higher needs to be installed on your PC. Additionally, it is 
recommended to install the package in a new environment.

The following commands are run in the command line.

1: Set up a new environment.
    
    python3 -m venv EmulsiPred_env

2: Enter (activate) the environment.

    source EmulsiPred_env/bin/activate

3a: Install EmulsiPred within the activated environment with pip.

    pip install EmulsiPred
    
3b: Install EmulsiPred by cloning from github and then installing with pip.

    git clone https://github.com/MarcatiliLab/EmulsiPred.git
    cd EmulsiPred
    pip install .
    cd ..
    

After either running 3a or 3b EmulsiPred is installed within the
environment you have created (in our case EmulsiPred_env).
---
#### Running EmulsiPred

After installation, EmulsiPred can be run from the terminal or
within a python script (with a little more flexibility from a
script).

EmulsiPred requires 2 inputs.
1) A fasta file containing the protein sequences to check for emulsifiers (termed sequence.fsa).
2) A NetSurfP file containing secondary structure information of the sequences in sequence.fsa (termed netsurfp.txt)  

Additionally, there are also 3 variable parameters. 
1) o: Output directory (default is current directory).
2) nr_seq: minimum number of sequences a peptide in the output is present in (default 1).
3) ls (lower score): minumum score a peptide in the output has (default 2).  

EmulsiPred can be run directly in the terminal with the following
command.

    python -m EmulsiPred -s path/to/sequence.fsa -n path/to/netsurfp.txt -o path/to/out_dir --nr_seq 1 --ls 2
    
Furthermore, it can be imported and run in a python script.

~~~~~~~~~~~~~~~~~~~~~
import EmulsiPred as ep

ep.EmulsiPred(sequences='path/to/sequence.fsa', netsurfp_results='path/to/netsurfp.txt', out_dir='path/to/out_dir', nr_seq=1, lower_score=2)
~~~~~~~~~~~~~~~~~~~~~