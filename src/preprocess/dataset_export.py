from Bio import SeqIO


def save_sequences(file_path, sequences):
    file_name = "{}/sequences.fa".format(file_path)
    SeqIO.write(sequences, file_name, "fasta")
