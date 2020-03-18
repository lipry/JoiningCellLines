from itertools import combinations

from src.config.config import config


def get_combinations_index(cell_lines, comb_length):
    indices = range(0, len(cell_lines))
    comb = []
    for r in comb_length:
        comb.extend(combinations(indices, r))

    return comb

