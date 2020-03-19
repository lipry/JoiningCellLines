import logging
from pandas.testing import assert_frame_equal
from src.preprocess.dataset_import import get_regions, get_categorical_labels, get_epigenetic_data
from src.preprocess.utils import append_without_duplicates, fill_missing, get_type


def preprocess_mode_exec(c):
    logging.basicConfig(format='[%(asctime)s] - %(levelname)s - %(message)s', level=logging.DEBUG)
    logging.debug("PREPROCESSING MODE")

    root_path = c['import_path']
    cell_lines = c['cell_lines']

    # Importing regions for enhancers and promoters
    enhancers_regions, promoters_regions = get_regions(root_path)

    # Importing and converting labels of enhancers and promoters and join them in a single dataframe
    full_sequences = get_categorical_labels("{}/{}".format(root_path, "fantom_1000"))

    # Importing epigenetic data
    logging.debug("Importing epigenetic data for: {}".format(", ".join(cell_lines)))
    for l in cell_lines:
        logging.debug("Importing {} data".format(l))

        df_epi_enanchers, df_epi_promoters = get_epigenetic_data("{}/{}".format(root_path, "fantom_1000"), l)

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

        # Check if the data are aligned dataframe are equals.
        assert len(full_sequences) == len(full_epi)
        assert_frame_equal(full_sequences[['chrom', 'chromStart', 'chromEnd']],
                           full_epi[['chrom', 'chromStart', 'chromEnd']])
