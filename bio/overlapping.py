import numpy as np
from .composition import string_composition

def overlapping_patterns(patterns: [str]) -> dict:
    """
    Input: A collection Patterns of k-mers.
    Output: The overlap graph Overlap(Patterns), in the form of an
            adjacency list.
    """
    adjacent = {}
    for i, pattern in enumerate(patterns):
        suffix = pattern[1:]
        for j, pattern_prime in enumerate(patterns):
            if i == j:
                continue

            if pattern_prime.startswith(suffix):
                if pattern not in adjacent:
                    adjacent[pattern] = []
                adjacent[pattern].append(pattern_prime)

    return adjacent


def main():
    #patterns = "ATGCG GCATG CATGC AGGCA GGCAT".split()
    with open("bio/data/dataset_198_9.txt") as f:
        patterns = f.read().strip().split()
        adjacent = overlapping_patterns(patterns)
        for key in sorted(adjacent):
            values = adjacent[key]
            for value in values:
                print("{0} -> {1}".format(key, value))


if __name__ == '__main__':
    main()


