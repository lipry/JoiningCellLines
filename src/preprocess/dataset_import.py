import pandas as pd

def import_labels(path):
    return pd.read_csv(path, sep="\t")


