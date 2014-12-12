import numpy as np
from .composition import string_composition
import random


def parse_graph(text: str) -> dict:
    output = text.strip().replace(" -> ", ":").split()
    graph = {}
    for x in output:
        key, values = x.split(":")
        if key not in graph:
            graph[key] = []
        graph[key] += sorted(values.split(","))
    return graph


def path_to_graph(path: list) -> dict:
    graph = {}
    for i, x in enumerate(path):
        if i < len(path) - 1 and x not in path:
            graph[x] = path[i + 1]
    return graph


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
    kmers = [text[i:i + k] for i in range(len(text) - k + 1)]
    return de_bruijn_graph(kmers)


def de_bruijn_graph(patterns: [str]) -> dict:
    """
    Input: A collection of k-mers Patterns.
    Output: The adjacency list of the de Bruijn graph DeBruijn(Patterns).
    """
    graph = {}
    for kmer in patterns:
        prefix_node, suffix_node = kmer[:-1], kmer[1:]
        if prefix_node not in graph:
            graph[prefix_node] = []
        graph[prefix_node].append(suffix_node)
    return graph


def eulerian_cycle(graph: dict) -> list:
    edges = [random.choice(list(graph.keys()))]
    marks = set()
    path = []
    while edges:
        node = edges[-1]
        has_unmarked_edge = False
        vertices = list(graph[node])
        random.shuffle(vertices)
        for vertex in vertices:
            mark = (node, vertex)
            if mark not in marks:
                marks.add(mark)
                edges.append(vertex)
                has_unmarked_edge = True
                break
        if not has_unmarked_edge:
            path.append(edges.pop())


    return list(reversed(path))


def main():
    graph = """
    0 -> 3
    1 -> 0
    2 -> 1,6
    3 -> 2
    4 -> 2
    5 -> 4
    6 -> 5,8
    7 -> 9
    8 -> 7
    9 -> 6
    """
    print("->".join(eulerian_cycle_sketch(parse_graph(graph))))
    #with open("bio/data/dataset_203_2.txt") as f:
    #    graph = f.read().strip()
    #text = "TAATGCCATGGGATGTT"
    #adjacent = path_graph(3, text)
    #with open("bio/data/dataset_200_7.txt") as f:
    #    #kmers = "GAGG CAGG GGGG GGGA CAGG AGGG GGAG".split()
    #    kmers = f.read().strip().split()
    #    adjacent = de_bruijn_graph(kmers)
    #    for key in sorted(adjacent):
    #        values = adjacent[key]
    #        print("{0} -> {1}".format(key, ",".join(sorted(values))))


if __name__ == '__main__':
    main()


