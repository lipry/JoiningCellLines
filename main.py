import logging
import os

from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

from utils.datasets_export import save_epigenetic, save_labels, save_sequences
from utils.datasets_import import import_cell_lines_data
from utils.datasets_join import multiple_intersect, merge
from utils.maps import Maps


def creating_seqIO(chr_id, sequence_string):
    return SeqRecord(sequence_string,
                     id=chr_id,
                     name="",
                     description="")


def building_final_datasets(ids, maps):
    final_seqIO = []
    final_epigenetic = []
    final_labels = []

    for chr_id in tqdm(ids):
        final_seqIO.append(creating_seqIO(chr_id, maps.get_sequence(chr_id)))
        final_epigenetic.append(maps.get_epigenetic_data(chr_id))
        final_labels.append(maps.get_label(chr_id))

    return final_seqIO, final_epigenetic, final_labels


cell_lines = ["GM12878", "HelaS3", "HepG2", "K562"]
files_path = "input_data"
output_path = "output_data"

# TODO: fix merge, for now intersect hard-coded
operation = 'i'
if operation not in ['m', 'i']:
    raise ValueError("Operation not valid")

if not os.path.exists(files_path):
    raise FileNotFoundError("Files path not found")

if not os.path.exists(output_path):
    logging.debug("{} not found, folder will be created")
    os.makedirs(output_path)

logging.basicConfig(format='[%(asctime)s] - %(levelname)s - %(message)s', level=logging.DEBUG)

logging.debug("Importing data from cell lines...")
X_sequences_list, X_epigenetic_list, y_list = import_cell_lines_data(files_path, cell_lines)

logging.debug("Building chr:interval id maps...")
d = Maps(len(cell_lines))
d.build_maps(X_sequences_list, X_epigenetic_list, y_list)

chr_ids = [[idx for idx, _ in line_data] for line_data in X_sequences_list]

if operation == 'i':
    logging.debug("Intersecting...")
    joined_ids = multiple_intersect(chr_ids)
else:
    logging.debug("Merging...")
    joined_ids = merge(chr_ids)

for cell_line, ids_list in zip(cell_lines, chr_ids):
    logging.debug("list {}: {}".format(cell_line, len(ids_list)))
logging.debug("Joined: {}".format(len(joined_ids)))

logging.debug("Bulding new intersected dataset...")
seqIO, epigenetic, labels = building_final_datasets(joined_ids, d)

logging.debug("Saving files...")
logging.debug("Saving sequences...")
file_suffix = "intersected" if operation == 'i' else "merged"
save_sequences(output_path, seqIO, file_suffix)

logging.debug("Saving epigenetic and labels data...")
for idx, cell_line in enumerate(cell_lines):
    epi = [e[idx] for e in epigenetic]
    lab = [e[idx] for e in labels]
    save_epigenetic(output_path, cell_line, epi, file_suffix)
    save_labels(output_path, cell_line, lab, file_suffix)

logging.debug("Finished.")
