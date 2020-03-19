import numpy as np

def append_without_duplicates(first_df, second_df):
    return first_df.append(second_df, ignore_index=True).drop_duplicates(
        subset=['chrom', 'chromStart', 'chromEnd'], keep=False)


def join_sequences_labels(enhancers, promoters):
    enhancers = enhancers.replace([1, 0], ["A-E", "I-E"])
    promoters = promoters.replace([1, 0], ["A-P", "I-P"])
    return append_without_duplicates(enhancers, promoters)


def fill_missing(df, metric="median"):
    for c in df.columns[4:]:
        if metric == "median":
            filled = df[c].median()
        df[c] = df[c].fillna(filled)
    return df


def get_type(c):
    d = {'chrom': np.object, 'chromStart': np.int, 'chromEnd': np.int, 'strand': np.object}
    try:
        return d[c]
    except KeyError:
        return np.float
