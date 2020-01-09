import numpy as np
from Bio import SeqIO


def save_sequences(file_path, sequences, suffix):
    file_name = "{}/sequences_{}.fa".format(file_path, suffix)
    SeqIO.write(sequences, file_name, "fasta")


def save_epigenetic(file_path, cell_lines, epigenetic_data, suffix):
    file_name = "{}/{}_epigenetic_{}.txt".format(file_path, cell_lines, suffix)
    np.savetxt(file_name, epigenetic_data, delimiter="\t", fmt="%.2f")


def save_labels(file_path, cell_lines, labels, suffix):
    file_name = "{}/{}_labels_{}.txt".format(file_path, cell_lines, suffix)
    np.savetxt(file_name, labels, fmt="%s")
