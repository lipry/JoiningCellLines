import os

import numpy as np
from Bio import SeqIO
from tqdm import tqdm

cell_lines = ["GM12878", "HelaS3", "HepG2", "K562"]


def check_cell_line(cell_line):
    if cell_line not in cell_lines:
        raise ValueError("Illegal cell line.")


def check_file_exist(file_path):
    if not os.path.exists(file_path):
        raise FileNotFoundError("file {} not found".format(file_path))


def import_epigenetic_data(files_path, cell_line):
    check_cell_line(cell_line)

    epigenetic_data_path = "{}/{}_200bp_Data.txt".format(files_path, cell_line)
    check_file_exist(epigenetic_data_path)

    epigenetic_data = np.loadtxt(epigenetic_data_path)

    return epigenetic_data


def import_sequences(files_path, cell_line):
    check_cell_line(cell_line)

    seqences_path = "{}/{}.fa".format(files_path, cell_line)
    check_file_exist(seqences_path)

    with open(seqences_path) as f:
        sequences = [(s.id, s.seq) for s in SeqIO.parse(f, 'fasta')]

    return sequences


def import_labels(files_path, cell_line):
    labels_path = "{}/{}_200bp_Classes.txt".format(files_path, cell_line)
    check_file_exist(labels_path)

    with open(labels_path, "r") as f:
        labels = np.array([line.strip() for line in f.readlines()])

    return labels


def import_cell_lines_data(files_path, cell_lines):
    X_sequences_list = []
    X_epigenetic_list = []
    y_list = []

    for line_name in tqdm(cell_lines):
        sequences = import_sequences(files_path, line_name)
        epigentic = import_epigenetic_data(files_path, line_name)
        labels = import_labels(files_path, line_name)

        X_sequences_list.append(sequences)
        X_epigenetic_list.append(epigentic)
        y_list.append(labels)

    return X_sequences_list, X_epigenetic_list, y_list
