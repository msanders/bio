import numpy as np
from .composition import (
    string_composition, genome_path_string, paired_composition, parse_patterns,
    string_spelled_by_gapped_patterns
)
from .numberpattern import generate_kmers
import networkx as nx
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


def build_di_graph(nodes: dict) -> nx.DiGraph:
    graph = nx.DiGraph()
    for node, values in nodes.items():
        for vertex in values:
            graph.add_edge(node, vertex, label="{0}->{1}".format(node, vertex))
    return graph


def build_dict_graph(nodes: nx.DiGraph) -> dict:
    graph = {}
    for left, right in nodes.edges():
        if left not in graph:
            graph[left] = []
        graph[left].append(right)
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


def paired_de_bruijn_graph(pairs: [(str, str)]) -> dict:
    graph = {}
    for pair in pairs:
        prefix_node = tuple(kmer[:-1] for kmer in pair)
        suffix_node = tuple(kmer[1:] for kmer in pair)
        if prefix_node not in graph:
            graph[prefix_node] = []
        graph[prefix_node].append(suffix_node)
    return graph


def eulerian_cycle(graph: dict, start=None) -> list:
    """
    Adopted from http://www.ms.uky.edu/~lee/ma515fa10/euler.pdf
    """
    edges = [start or random.choice(list(graph.keys()))]
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


def eulerian_cycle2(graph: dict) -> list:
    graph = build_di_graph(graph)
    return [
        x[0] for x in nx.eulerian_circuit(graph, random.choice(graph.nodes()))
    ]


def find_degrees(graph: nx.DiGraph) -> (dict, dict):
    in_degrees = {}
    out_degrees = {}
    for left, values in graph.items():
        out_degrees[right] = len(values)
        for right in values:
            in_degrees[right] = in_degrees.get(right, 0) + 1
    return in_degrees, out_degrees


def in_out_balance(graph: nx.DiGraph, node: object):
    return graph.in_degree(node) - graph.out_degree(node)


def balance_graph(graph: nx.DiGraph) -> (object, object):
    head = None
    tail = None
    marked_nodes = set()
    for node in graph:
        balance = in_out_balance(graph, node)
        if abs(balance) > 1:
            raise ValueError("Cannot balance graph: {0}".format(graph))
        elif balance == 1:
            tail = node
        elif balance == -1:
            head = node

    edge = (tail, head)
    if head is not None:
        graph.add_edge(*edge)
    return edge


def euler_path(graph: nx.DiGraph, func=genome_path_string) -> list:
    edge = balance_graph(graph)
    if not nx.is_eulerian(graph):
        raise ValueError("Not Eulerian: {0}".format(graph))

    circuit = list(nx.eulerian_circuit(graph, edge[1]))
    #print("asdf {0}".format(circuit))
    #return [ func(x) for x in circuit ]
    return [x[0] for x in circuit] + [ circuit[0][0] ]


def euler_path2(graph: dict) -> list:
    nx_graph = build_di_graph(graph)
    edge = balance_graph(nx_graph)
    if not nx.is_eulerian(nx_graph):
        raise ValueError("Not Eulerian: {0}".format(graph))

    return eulerian_cycle(build_dict_graph(nx_graph), edge[1])


def reconstruct_string(patterns: [str], k: int) -> str:
    """
    Input: An integer k followed by a list of k-mers Patterns.
    Output: A string Text with k-mer composition equal to Patterns.
    """
    graph = de_bruijn_graph(patterns)
    graph = build_di_graph(graph)
    return genome_path_string(euler_path(graph)[:-1])


def binary_strings(k: int) -> [str]:
    return generate_kmers("01", k)


def is_universal(text: str, k) -> bool:
    kmers = list(binary_strings(k))
    rotated_text = text + text[:k - 1]
    count = 0
    for kmer in kmers:
        if rotated_text.count(kmer) != 1:
            return False
        count += 1
    return count == len(kmers)


def universal_circular_string(k: int) -> str:
    kmers = list(binary_strings(k))
    graph = de_bruijn_graph(kmers)
    path = euler_path(build_di_graph(graph))
    return genome_path_string(path[:-(k - 1)])

def brute_readpair_reconstruction(pairs: [(str, str)], d: int) -> str:
    for _ in range(5000):
        path = euler_path(build_di_graph(paired_de_bruijn_graph(pairs)))
        text = string_spelled_by_gapped_patterns(path, d)
        if text is not None:
            break
    return text


def main():
    #print(is_universal(3, "0001110100"))
    print(is_universal(4, "0000110010111101"))
    #print(is_universal(4, "0001111011001010000"))
    k = 4
    t = universal_circular_string(k)
    print(t)
    print(is_universal(k, t))
    #with open("bio/data/dataset_203_6.txt") as f:
    #    kmers = f.read().strip().split()
    #    print(reconstruct_string(kmers, 25))
    #print("->".join(eulerian_cycle_sketch(parse_graph(graph))))
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


