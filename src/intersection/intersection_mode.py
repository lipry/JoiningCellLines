import logging
import os
from operator import itemgetter

from src.intersection.datasets_export import save_epigenetic, save_labels, save_sequences, building_final_datasets
from src.intersection.datasets_import import import_cell_lines_data
from src.intersection.datasets_join import multiple_intersect, merge
from src.intersection.maps import Maps
from src.intersection.utils import get_combinations_index


def intersection_mode_exec(c):
    logging.basicConfig(format='[%(asctime)s] - %(levelname)s - %(message)s', level=logging.DEBUG)
    logging.debug("INTERSECTION MODE")

    cell_lines = c['cell_lines']
    files_path = c['import_path']
    output_path = c['export_path']

    # TODO: fix merge, for now intersect hard-coded
    # TODO: refactor, removing merge.
    operation = 'i'
    if operation not in ['m', 'i']:
        raise ValueError("Operation not valid")

    if not os.path.exists(files_path):
        raise FileNotFoundError("Files path not found: {}".format(files_path))

    if not os.path.exists(output_path):
        logging.debug("{} not found, folder will be created")
        os.makedirs(output_path)



    logging.debug("Importing data from cell lines...")
    X_sequences_list, X_epigenetic_list, y_list = import_cell_lines_data(files_path, cell_lines)

    comb = get_combinations_index(cell_lines, c['combinations_length'])

    for lines_idx in comb:
        # Filter combinations
        cell_line_comb = list(itemgetter(*lines_idx)(cell_lines))
        X_sequences_list_comb = list(itemgetter(*lines_idx)(X_sequences_list))
        X_epigenetic_list_comb = list(itemgetter(*lines_idx)(X_epigenetic_list))
        y_list_comb = list(itemgetter(*lines_idx)(y_list))

        logging.debug("NEW COMBINATION: {}".format(cell_line_comb))

        logging.debug("Building chr:interval id maps...")
        d = Maps(len(cell_line_comb))
        d.build_maps(X_sequences_list_comb, X_epigenetic_list_comb, y_list_comb)

        chr_ids = [[idx for idx, _ in line_data] for line_data in X_sequences_list_comb]

        if operation == 'i':
            logging.debug("Intersecting...")
            joined_ids = multiple_intersect(chr_ids)
        else:
            logging.debug("Merging...")
            joined_ids = merge(chr_ids)

        for cell_line, ids_list in zip(cell_line_comb, chr_ids):
            logging.debug("list {}: {}".format(cell_line, len(ids_list)))
        logging.debug("Joined: {}".format(len(joined_ids)))

        logging.debug("Bulding new intersected dataset...")
        seqIO, epigenetic, labels, diff_label_counter = building_final_datasets(joined_ids, d)
        logging.debug("Numeber of elements with different labels in intersections: {}".format(diff_label_counter))

        if c['save_results']:
            logging.debug("Saving files...")
            logging.debug("Saving sequences...")
            file_suffix = "intersected" if operation == 'i' else "merged"
            save_sequences(output_path, seqIO, file_suffix)

            logging.debug("Saving epigenetic and labels data...")
            for idx, cell_line in enumerate(cell_line_comb):
                epi = [e[idx] for e in epigenetic]
                lab = [e[idx] for e in labels]
                save_epigenetic(output_path, cell_line, epi, file_suffix)
                save_labels(output_path, cell_line, lab, file_suffix)

    logging.debug("Finished.")
