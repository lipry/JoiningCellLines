import heapq

import itertools


def intersect(chr1, chr2):
    results = []
    p1 = 0
    p2 = 0
    while p1 < len(chr1) and p2 < len(chr2):
        if chr1[p1] == chr2[p2]:
            results.append(chr1[p1])
            p1 += 1
            p2 += 1
        else:
            if chr1[p1] < chr2[p2]:
                p1 += 1
            else:
                p2 += 1
    return results


def sort_cell_lines(chr_list):
    return sorted([sorted(chromosom) for chromosom in chr_list], key=lambda x: len(x))


def multiple_intersect(chr_list):
    cell_lines = sort_cell_lines(chr_list)

    results = cell_lines[0]
    for line in cell_lines[1:]:
        results = intersect(results, line)
    return results


def merge(chr_list):
    return sorted(list(set(itertools.chain(*chr_list))))
