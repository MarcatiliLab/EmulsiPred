# EmulsiPred
Prediction of Emulsifying Peptides


###Prerequisites and installation

(This part will be written soon)



python-3.6 or higher needs to be installed on your PC.

move EmulsiPred folder to preferred destination.

Install the package by running setup.py. This is done with the
following command.

    pip install .
    


Predict emulsifying peptides with the following.

    python EmulsiPred.py -s ../Examples/examples.fsa -n ../Examples/example_netsurfp.txt -v ../NormalizationValues/ -o ../Results/ --n_seq 1 --ls 3