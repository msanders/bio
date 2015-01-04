import numpy as np
import random
import re
from .graph import build_di_graph
from copy import deepcopy


def parse_graph(text: str, weighted: bool = False) -> dict:
    output = text.strip().replace(" -> ", "->").split()
    graph = {}
    for x in output:
        key, values = x.split("->")
        if weighted:
            key = int(key)
        if key not in graph:
            graph[key] = []
        if weighted:
            graph[key] = sorted(graph[key] + [tuple(
                int(x) for x in values.split(":"))]
            )
        else:
            graph[key] = sorted(graph[key] + values.split(","))
    return graph


def manhatten_tourist(n: int, m: int, down: np.ndarray, right: np.ndarray):
    s = np.zeros([n + 1, m + 1], dtype='i')
    s[0][0] = 0
    for i in range(1, n + 1):
        s[i][0] = s[i - 1][0] + down[i - 1][0]

    for j in range(1, m + 1):
        s[0][j] = s[0][j - 1] + right[0][j - 1]

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s[i][j] = max(s[i - 1][j] + down[i-1][j], s[i][j - 1] + right[i][j-1])
    return s[n][m]


def parse_matrix(text: str) -> np.ndarray:
    return np.asarray([
        [int(x) for x in line.split()] for line in text.strip().splitlines()
    ], dtype='i')


def parse_tourist_input(text: str) -> (int, int, np.ndarray, np.ndarray):
    lines = text.strip().splitlines()
    n, m = lines[0].split()
    n, m = int(n), int(m)
    down = parse_matrix("\n".join(lines[1:][:n]))
    right = parse_matrix("\n".join(lines[2 + n:]))
    return n, m, down, right


def lcs_path(v: str, w: str) -> np.ndarray:
    s = np.zeros((len(v) + 1, len(w) + 1), dtype='i')
    for i in range(len(v)):
        for j in range(len(w)):
            if v[i] == w[j]:
                s[i + 1][j + 1] = s[i][j] + 1
            else:
                s[i + 1][j + 1] = max(s[i + 1][j], s[i][j + 1])
    return s


def lcs_scored_backtrack(matrix: dict, v: str, w: str, o: int, u: int, min_cost: int) -> (
    np.ndarray, np.ndarray
):
    s = np.zeros((len(v) + 1, len(w) + 1), dtype='i')
    backtrack = np.zeros((len(v) + 1, len(w) + 1), dtype='str')

    for i in range(1, len(v) + 1):
        s[i][0] = -i * min_cost
    for j in range(1, len(w) + 1):
        s[0][j] = -j * min_cost

    for i in range(1, len(v) + 1):
        for j in range(1, len(w) + 1):
            s[i][j], backtrack[i][j] = max(
                (s[i - 1][j] - o, "↓"),
                (s[i][j - 1] - o, "→"),
                (s[i - 1][j - 1] + matrix[v[i - 1]][w[j - 1]], "↘"),
                (-i * j * min_cost, "_"),
                key=lambda x: x[0]
            )
    return s, backtrack


def lcs_backtrack(s: np.ndarray) -> np.ndarray:
    backtrack = np.zeros([len(s), len(s[0])], dtype='str')

    for i in range(1, len(s)):
        for j in range(1, len(s[0])):
            if s[i][j] == s[i - 1][j]:
                backtrack[i][j] = "↓"
            elif s[i][j] == s[i][j - 1]:
                backtrack[i][j] = "→"
            else:
                backtrack[i][j] = "↘"
    return backtrack


def output_lcs(backtrack: np.ndarray, v: str, w: str):
    lcs = []
    i, j = len(v), len(w)
    while i > 0 and j > 0:
        if backtrack[i][j] == "↓":
            i -= 1
        elif backtrack[i][j] == "→":
            j -= 1
        elif backtrack[i][j] == "↘":
            lcs.append(v[i - 1])
            i -= 1
            j -= 1
    return "".join(reversed(lcs))


def outgoing_edges(graph: dict) -> set:
    return set(x for values in graph.values() for x in values)


def topological_ordering(graph: dict) -> list:
    graph = deepcopy(graph)
    order = []
    values = outgoing_edges(graph)
    candidates = set(graph.keys()) - values

    while candidates:
        node = random.sample(candidates, 1)[0]
        order.append(node)
        candidates.remove(node)

        outgoing_nodes = graph.get(node, []).copy()
        graph[node] = []
        values = outgoing_edges(graph)

        for node in outgoing_nodes:
            if node not in values:
                candidates.add(node)

    if outgoing_edges(graph):
        raise ValueError("Input graph is not a DAG")
    return order


def unweighted_graph(graph: dict) -> dict:
    unweighted = {}
    for key in graph:
        unweighted[key] = [x[0] for x in graph[key]]
    return unweighted


def incoming_edges(graph: dict, node: object) -> list:
    return [key for key in graph if node in graph[key]]


def edge_with_weight(graph: dict, in_node: object, out_node: object) -> int:
    return [x for x in graph[in_node] if x[0] == out_node][0]


def longest_dag_path(graph: dict, source: int, sink: int) -> (int, list):
    unweighted = unweighted_graph(graph)
    s = {node: float("-inf")
         for node in set(graph) | outgoing_edges(unweighted)}
    s[source] = 0

    top_order = topological_ordering(unweighted)
    top_order = top_order[top_order.index(source) + 1:top_order.index(sink) + 1]

    backtrack = {node: None for node in top_order}
    for node in top_order:
        for key in incoming_edges(unweighted, node):
            node, weight = edge_with_weight(graph, key, node)
            weight += s[key]
            if weight > s[node]:
                s[node], backtrack[node] = weight, key

    # Backtrack to get the longest path.
    path = [sink]
    while path[0] != source:
        path.insert(0, backtrack[path[0]])

    return s[sink], path


def parse_matrix(text: str) -> dict:
    matrix = {}
    lines = text.strip().splitlines()
    keys = [x for x in lines[0].split()]

    for key in keys:
        matrix[key] = {}

    for xkey, line in zip(keys, lines[1:]):
        values = line.split()[1:]
        for ykey, value in zip(keys, values):
            matrix[xkey][ykey] = int(value)

    return matrix


def blosum62():
    with open("bio/data/BLOSUM62.txt") as f:
        return parse_matrix(f.read())


def pam250():
    with open("bio/data/PAM250.txt") as f:
        return parse_matrix(f.read())


def compute_alignment(matrix: dict, v: str,
                      w: str, o: int, u: int,
                      min_cost: int,
                      backtrack_sink: object,
                      backtrack_while: object) -> (int, str, str):
    s, backtrack = lcs_scored_backtrack(matrix, v, w, o, u, min_cost)

    v_aligned = ""
    w_aligned = ""

    i, j = backtrack_sink(s)
    max_score = s[i][j]
    while backtrack_while(backtrack, i, j):
        direction = backtrack[i][j]
        if direction == "↓":
            v_aligned += v[i - 1]
            w_aligned += "-"
            i -= 1
        elif direction == "→":
            v_aligned += "-"
            w_aligned += w[j - 1]
            j -= 1
        elif direction == "↘":
            v_aligned += v[i - 1]
            w_aligned += w[j - 1]
            i -= 1
            j -= 1

    # Prepend necessary indels to get back to (0, 0).
    if "_" not in backtrack.flatten():
        for repeat in range(i):
            w_aligned += "-"
            v_aligned += v[i - 1]

        for repeat in range(j):
            v_aligned += "-"
            w_aligned += w[j - 1]

    return (
        max_score,
        "".join(reversed(v_aligned)),
        "".join(reversed(w_aligned))
    )


def global_alignment_problem(matrix: dict, v: str, w: str, o: int, u: int) -> (
    int, str, str
):
    return compute_alignment(
        matrix, v, w, o, u,
        min_cost=o,
        backtrack_sink=lambda _: (len(v), len(w)),
        backtrack_while=lambda _, i, j: i > 0 and j > 0
    )


def local_alignment_problem(matrix: dict, v: str, w: str, o: int, u: int) -> (
    int, str, str
):
    return compute_alignment(
        matrix, v, w, o, u,
        min_cost=0,
        backtrack_sink=lambda s: np.unravel_index(s.argmax(), s.shape),
        backtrack_while=lambda backtrack, i, j: i > 0 and j > 0 and backtrack[i][j] != "_"
    )


def levenshtein_distance(v: str, w: str) -> int:
    """
    From http://hetland.org/coding/python/levenshtein.py
    """
    n, m = len(v), len(w)
    if n > m:
        # Make sure n <= m, to use O(min(n,m)) space
        v,w = w,v
        n,m = m,n

    current = range(n + 1)
    for i in range(1, m + 1):
        previous, current = current, [i] + [0] * n
        for j in range(1, n + 1):
            add, delete = previous[j] + 1, current[j - 1] + 1
            change = previous[j - 1]
            if v[j - 1] != w[i - 1]:
                change = change + 1
            current[j] = min(add, delete, change)

    return current[n]


def main():
    print(levenshtein_distance("PLEASANTLY", "MEANLY"))
    print(global_alignment_problem(blosum62(), "PLEASANTLY", "MEANLY", u=0, o=5))
    print(local_alignment_problem(pam250(), "MEANLY", "PENALTY", u=0, o=5))
    return
    #text = """
    #0 -> 1,11,12,14,15,16,17,18,2,3,6,7,8
    #1 -> 11,14,17,18,19,2,7,8
    #10 -> 12,13,14,15,18,19
    #11 -> 13,14,17
    #12 -> 15,16,19
    #13 -> 18
    #14 -> 15,16,17,19
    #15 -> 17,19
    #16 -> 17
    #17 -> 19
    #2 -> 14,15,19,3,9
    #3 -> 14,15,17,4,5,7,8,9
    #4 -> 11,13,14,15,18,6,8
    #5 -> 10,12,6,8,9
    #6 -> 10,12,13,15,17,18,19,7
    #7 -> 10,12,13,14,16,17,19,9
    #8 -> 10,13,14,15,16,18,9
    #9 -> 12,14,17
    #"""

    #graph = parse_graph(text)
    #print(", ".join(topological_ordering(graph)))
    #return
    text = """
    0->1:7
    0->2:4
    2->3:2
    1->4:1
    3->4:3
    """
    graph = parse_graph(text, weighted=True)
    print(graph)
    print(longest_dag_path(graph, 0, 4))
    return
    text = """AACCTTGG
    ACACTGTGA"""
    v, w = text.split()
    s = lcs_path(v, w)
    backtrack = lcs_backtrack(s)
    print(backtrack)
    lcs = output_lcs(backtrack, s, v, w)
    print(lcs)
    return

    text = """
    6 15
    2 0 4 2 0 1 4 4 1 1 2 4 2 3 3 3
    0 3 2 2 0 1 4 1 1 3 1 3 3 4 3 3
    4 1 0 3 1 2 1 0 3 1 2 3 2 0 4 4
    4 1 1 4 2 1 3 4 0 3 4 1 3 3 4 0
    4 1 4 4 0 3 1 0 3 0 3 3 3 0 4 4
    1 0 3 0 0 1 4 4 0 2 2 3 1 2 0 2
    -
    1 3 1 3 4 0 1 2 0 4 1 0 0 2 0
    0 3 1 1 0 2 4 4 4 1 1 2 1 1 3
    3 3 1 1 2 0 4 0 4 1 4 3 2 4 2
    4 0 2 1 0 0 0 1 1 2 4 4 4 3 3
    2 2 3 3 4 2 1 3 2 3 4 3 1 0 1
    0 2 3 3 0 3 4 0 3 4 3 3 2 0 0
    3 4 4 1 0 1 3 0 2 3 2 2 2 0 3
    """
    n, m, down, right = parse_tourist_input(text)
    print(manhatten_tourist(n, m, down, right))


if __name__ == '__main__':
    main()
