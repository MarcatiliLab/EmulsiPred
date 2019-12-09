from Bio import ExPASy
from Bio import SwissProt
import argparse
import numpy as np
from openpyxl import load_workbook
from openpyxl.utils import column_index_from_string as colindex
from openpyxl.utils import coordinate_to_tuple as cototu
import urllib.error

## If you want to find sequences for the accession numbers in fasta
# python Acc_to_seq.py -e /home/projects/vaccine/people/s132421/Arbejde/Data/patatin_and_kunitz.xlsx -o /home/projects/vaccine/people/s132421/Arbejde/Data

parser = argparse.ArgumentParser()
parser.add_argument('-e', dest='accession_list', default='', help='Excel file with accession numbers')
parser.add_argument('-o', dest='out_dir', default='', help='Path to directory output will be saved in')

args = parser.parse_args()
excelfile = args.accession_list
out_dir = args.out_dir


def get_accession(xfi):

    acclst = []  # list to hold accession numbers
    if xfi != '':
        wbt = load_workbook(filename=xfi, read_only=True)
        wst = wbt.active
        for row in wst.rows:
            for cell in row:
                if cell.value:
                    if cototu(cell.coordinate)[1] == 1:
                        a = str(cell.value)
                        acclst.append(a.strip())
                    else:
                        break
    return(acclst)

def seq_from_accession(accession_nr):
    try:
        handle = ExPASy.get_sprot_raw(accession_nr)
    except urllib.error.HTTPError:
        print("problem with:", accession_nr)
        return

    record = SwissProt.read(handle)
    return(record.sequence)

def seq_from_accession_list(accession_excel_file):

    # First take out accession numbers from the excel file and put in a put
    acclst = get_accession(accession_excel_file)
    seqs = {}

    for accession in acclst:

        # Second retrieve the sequence from each accession number
        seq = seq_from_accession(accession)
        seqs[accession] = seq

    return(seqs)

def fasta_files(excelfile):

    sequences = seq_from_accession_list(excelfile)
    print(sequences)

    with open(out_dir+'/PeptideSequences.fsa','w') as out_file:
        for acces in sequences:
            out_file.write('>' + acces + '\n')
            out_file.write(sequences[acces] + '\n')

    return

fasta_files(excelfile)


