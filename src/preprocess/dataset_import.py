import gzip

import pandas as pd

from src.preprocess.utils import join_sequences_labels


def _check_type(t):
    if t not in ['enhancers', 'promoters']:
        raise ValueError("Type must be enhancers or promoters, {} found".format(t))


def _import_region(root_path, t="enhancers"):
    _check_type(t)

    return pd.read_csv("{}/{}_regions.bed".format(root_path, t),
                       sep="\t",
                       names=['chrom', 'chromStart', 'chromEnd', 'full_code', 'unk', 'unk1'])


def get_regions(root_path):
    enhancers = _import_region(root_path, t="enhancers")
    promoters = _import_region(root_path, t="promoters")

    return enhancers, promoters


def _import_binary_labels(root_path, t="enhancers"):
    _check_type(t)
    return pd.read_csv("{}/{}.bed".format(root_path, t), sep="\t")


def get_binary_labels(root_path):
    enhancers = _import_binary_labels(root_path, t="enhancers")
    promoters = _import_binary_labels(root_path, t="promoters")

    return enhancers, promoters


def get_categorical_labels(root_path):
    en_bin, prm_bin = get_binary_labels(root_path)
    return join_sequences_labels(en_bin, prm_bin)


def import_epigenetic_data(root_path, cell_line, t="enhancers"):
    _check_type(t)

    with gzip.open('{}/{}/{}.csv.gz'.format(root_path, t, cell_line)) as f:
        df_epi = pd.read_csv(f, low_memory=False)
    return df_epi[1:]


def get_epigenetic_data(root_path, cell_line):
    enhancers_epi = import_epigenetic_data(root_path, cell_line, t="enhancers")
    promoters_epi = import_epigenetic_data(root_path, cell_line, t="promoters")

    return enhancers_epi, promoters_epi
