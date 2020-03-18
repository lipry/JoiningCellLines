import gzip
import logging
import numpy as np

import pandas as pd
from pandas.util.testing import assert_frame_equal


def get_type(c):
    d = {'chrom': np.object, 'chromStart': np.int, 'chromEnd': np.int, 'strand': np.object}
    try:
        return d[c]
    except KeyError:
        return np.float

def preprocess_mode_exec(c):
    logging.basicConfig(format='[%(asctime)s] - %(levelname)s - %(message)s', level=logging.DEBUG)
    logging.debug("PREPROCESSING MODE")

    root_path = c['import_path']
    cell_lines = c['cell_lines']

    # Importing regions for enhancers and promoters
    enhancers_regions = pd.read_csv(root_path + "/enhancers_regions.bed",
                                    sep="\t",
                                    names=['chrom', 'chromStart', 'chromEnd', 'full_code', 'unk', 'unk1'])
    promoters_regions = pd.read_csv(root_path + "/promoters_regions.bed",
                                    sep="\t",
                                    names=['chrom', 'chromStart', 'chromEnd', 'full_code', 'unk', 'unk1'])

    # Importing and converting labels of enhancers and promoters
    enhancers = pd.read_csv(root_path + "/fantom_1000/enhancers.bed", sep="\t")
    promoters = pd.read_csv(root_path + "/fantom_1000/promoters.bed", sep="\t")
    enhancers.replace([1, 0], ["A-E", "I-E"], inplace=True)
    promoters.replace([1, 0], ["A-P", "I-P"], inplace=True)
    full_sequences = enhancers.append(promoters, ignore_index=True)\
                        .drop_duplicates(subset =['chrom', 'chromStart', 'chromEnd'],
                                         keep = False)

    #unire enhancers e promoters

    # Importing epigenetic data
    logging.debug("Importing epigenetic data for: {}".format(", ".join(cell_lines)))
    for l in cell_lines:
        logging.debug("Importing {} data".format(l))

        with gzip.open(root_path + '/fantom_1000/enhancers/{}.csv.gz'.format(l)) as f:
            df_epi_enanchers = pd.read_csv(f,
                                           low_memory=False)
        df_epi_enanchers = df_epi_enanchers[1:]


        with gzip.open(root_path + '/fantom_1000/promoters/{}.csv.gz'.format(l)) as f:
            df_epi_promoters = pd.read_csv(f,
                                           low_memory=False)
        df_epi_promoters = df_epi_promoters[1:]

        # building type dictionary
        converting_dictionary = {c: get_type(c) for c in df_epi_promoters.columns}
        df_epi_promoters = df_epi_promoters.astype(converting_dictionary)
        df_epi_enanchers = df_epi_enanchers.astype(converting_dictionary)

        assert len(df_epi_promoters.columns) == len(df_epi_enanchers.columns)
        logging.debug("number features for {}: {}".format(l, len(df_epi_promoters.columns)-4))
        logging.debug("Number of missing values in promoters: {}".format(df_epi_promoters.isna().sum().sum()))
        logging.debug("Number of missing values in enhancers: {}".format(df_epi_enanchers.isna().sum().sum()))

        for c in df_epi_promoters.columns[4:]:
            df_epi_promoters[c].fillna(df_epi_promoters[c].median(), inplace=True)

        for c in df_epi_enanchers.columns[4:]:
            df_epi_enanchers[c].fillna(df_epi_enanchers[c].median(), inplace=True)

        assert len(enhancers_regions) == len(df_epi_enanchers)
        logging.debug("Enhancers - regions: {}, epigenetics: {}".format(len(enhancers_regions),
                                                                        len(df_epi_enanchers)))

        assert len(promoters_regions) == len(df_epi_promoters)
        logging.debug("Promoters - regions: {}, epigenetics: {}".format(len(promoters_regions),
                                                                        len(df_epi_promoters)))

        full_epi = df_epi_enanchers.append(df_epi_promoters, ignore_index=True)\
                        .drop_duplicates(subset =['chrom', 'chromStart', 'chromEnd'],
                                         keep = False)

        # Check if the data are aligned dataframe are equals.
        assert len(full_sequences) == len(full_epi)
        assert_frame_equal(full_sequences[['chrom', 'chromStart', 'chromEnd']],
                           full_epi[['chrom', 'chromStart', 'chromEnd']])