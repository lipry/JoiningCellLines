import numpy as np
from tqdm import tqdm
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def save_sequences(file_path, sequences, suffix):
    file_name = "{}/sequences_{}.fa".format(file_path, suffix)
    SeqIO.write(sequences, file_name, "fasta")


def save_epigenetic(file_path, cell_lines, epigenetic_data, suffix):
    file_name = "{}/{}_epigenetic_{}.txt".format(file_path, cell_lines, suffix)
    np.savetxt(file_name, epigenetic_data, delimiter="\t", fmt="%.2f")


def save_labels(file_path, cell_lines, labels, suffix):
    file_name = "{}/{}_labels_{}.txt".format(file_path, cell_lines, suffix)
    np.savetxt(file_name, labels, fmt="%s")


def creating_seqIO(chr_id, sequence_string):
    return SeqRecord(sequence_string,
                     id=chr_id,
                     name="",
                     description="")


def building_final_datasets(ids, maps):
    final_seqIO = []
    final_epigenetic = []
    final_labels = []
    different_labels_counts = 0

    for chr_id in tqdm(ids):
        final_seqIO.append(creating_seqIO(chr_id, maps.get_sequence(chr_id)))
        final_epigenetic.append(maps.get_epigenetic_data(chr_id))
        final_labels.append(maps.get_label(chr_id))
        if not all(x==maps.get_label(chr_id)[0] for x in maps.get_label(chr_id)):
            different_labels_counts+=1

    return final_seqIO, final_epigenetic, final_labels, different_labels_counts
