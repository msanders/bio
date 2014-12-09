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


def de_bruijn_path(k: int, text: str) -> dict:
    graph = {}
    for i in range(len(text) - k + 1):
        kmer = text[i:i + k]
        prefix_node, suffix_node = kmer[:-1], kmer[1:]
        if prefix_node not in graph:
            graph[prefix_node] = []
        graph[prefix_node].append(suffix_node)
    return graph


def main():
    #text = "TAATGCCATGGGATGTT"
    #adjacent = path_graph(3, text)
    k = 4
    text = f.read().strip()
    adjacent = de_bruijn_path(k, text)
    for key in sorted(adjacent):
        values = adjacent[key]
        print("{0} -> {1}".format(key, ",".join(values)))


if __name__ == '__main__':
    main()


