import logging
import os

import numpy as np
from Bio.Seq import Seq
from pandas.testing import assert_frame_equal
from ucsc_genomes_downloader import Genome

from src.intersection.datasets_export import creating_seqIO
from src.preprocess.dataset_export import save_sequences
from src.preprocess.dataset_import import get_regions, get_categorical_labels, get_epigenetic_data, get_full_path
from src.preprocess.utils import append_without_duplicates, fill_missing, get_type


def preprocess_mode_exec(c):
    logging.basicConfig(format='[%(asctime)s] - %(levelname)s - %(message)s', level=logging.DEBUG)
    logging.debug("PREPROCESSING MODE")

    root_path = c['import_path']
    saving_path = c['export_path']
    cell_lines = c['cell_lines']
    window_size = c['window_size']
    dataset_type = c['dataset']


    if not os.path.exists(root_path):
        raise FileNotFoundError("Files path not found: {}".format(root_path))

    if not os.path.exists(saving_path):
        logging.debug("{} not found, folder will be created")
        os.makedirs(saving_path)

    label_epi_path = get_full_path(root_path, window_size, dataset_type)

    # Importing regions for enhancers and promoters
    enhancers_regions, promoters_regions = get_regions(root_path)

    # Importing and converting labels of enhancers and promoters and join them in a single dataframe
    full_sequences = get_categorical_labels(label_epi_path)
    logging.debug("Saving the sequences bed file in {}".format(saving_path))
    full_sequences.to_csv("{}/sequences.bed".format(saving_path),
                          sep="\t",
                          columns=['chrom', 'chromStart', 'chromEnd'],
                          header=False,
                          index=False)

    logging.debug("Downloading the hg19 genome")
    chroms = [k for k, _ in full_sequences.groupby(['chrom'])]
    hg19 = Genome(assembly="hg19", chromosomes=chroms)
    logging.debug("Downloading the hg19 genome")
    sequences = hg19.bed_to_sequence(full_sequences)

    logging.debug(sequences.loc[(sequences['chrom'] == "chr1") & (sequences['chromStart'] == 839287) & (sequences['chromEnd'] == 840287)])

    logging.debug("Saving sequences to file...")
    seqIO_seq = [creating_seqIO("{}:{}-{}".format(row['chrom'], row['chromStart'], row['chromEnd']), Seq(row['sequence'].upper()))
                 for _, row in sequences.iterrows()]
    save_sequences(saving_path, seqIO_seq)


    # Importing epigenetic data
    logging.debug("Importing epigenetic data for: {}".format(", ".join(cell_lines)))
    logging.debug("-------------------------------------------------------------")
    for l in cell_lines:
        logging.debug("Importing {} data".format(l))

        df_epi_enanchers, df_epi_promoters = get_epigenetic_data(label_epi_path, l)

        # building type dictionary
        converting_dictionary = {c: get_type(c) for c in df_epi_promoters.columns}
        df_epi_enanchers = df_epi_enanchers.astype(converting_dictionary)
        df_epi_promoters = df_epi_promoters.astype(converting_dictionary)

        assert len(df_epi_promoters.columns) == len(df_epi_enanchers.columns)
        logging.debug("number features for {}: {}".format(l, len(df_epi_promoters.columns)-4))
        logging.debug("Number of missing values in enhancers: {}".format(df_epi_enanchers.isna().sum().sum()))
        logging.debug("Number of missing values in promoters: {}".format(df_epi_promoters.isna().sum().sum()))

        df_epi_enanchers = fill_missing(df_epi_enanchers, metric="median")
        df_epi_promoters = fill_missing(df_epi_promoters, metric="median")

        assert len(enhancers_regions) == len(df_epi_enanchers)
        logging.debug("Enhancers - regions: {}, epigenetics: {}".format(len(enhancers_regions),
                                                                        len(df_epi_enanchers)))

        assert len(promoters_regions) == len(df_epi_promoters)
        logging.debug("Promoters - regions: {}, epigenetics: {}".format(len(promoters_regions),
                                                                        len(df_epi_promoters)))

        full_epi = append_without_duplicates(df_epi_enanchers, df_epi_promoters)

        # Check if the data are aligned dataframe are equals before save.
        assert len(full_sequences) == len(full_epi)
        assert_frame_equal(full_sequences[['chrom', 'chromStart', 'chromEnd']],
                           full_epi[['chrom', 'chromStart', 'chromEnd']])
        logging.debug("Number of total sequences: {}".format(len(full_sequences)))

        logging.debug("Saving results in {}".format(saving_path))
        np.savetxt("{}/{}_epigenetic.txt".format(saving_path, l), full_epi.iloc[: , 4:].values, fmt='%f')
        np.savetxt("{}/{}_labels.txt".format(saving_path, l), full_sequences[l].values, fmt='%s')

        logging.debug("-------------------------------------------------------------")
