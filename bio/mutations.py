import pyximport; pyximport.install()
from .chamming import hamming_distance
from .numberpattern import all_nucleotides
from .reverse_complement import complement
import itertools

NUCLEOTIDES = "ATCG"


def mutate_complements(pattern):
    mutations = [pattern]
    for i in range(len(pattern)):
        mutation = bytearray(pattern, "utf8")
        mutation[i] = ord(complement(pattern[i]))
        mutations.append("".join(map(chr, mutation)))
    return mutations


def immediate_neighbors(pattern):
    mutations = {pattern}
    for i in range(len(pattern)):
        for nucleotide in NUCLEOTIDES:
            if nucleotide != pattern[i]:
                mutation = bytearray(pattern, "utf8")
                mutation[i] = ord(nucleotide)
                mutations.add("".join(map(chr, mutation)))
    return mutations


def neighbors_recursive(pattern, d):
    if d == 0:
        return [pattern]
    elif len(pattern) == 1:
        return {"A", "C", "G", "T"}

    neighborhood = set()
    suffix_neighbors = neighbors(pattern[1:], d)
    for text in suffix_neighbors:
        if hamming_distance(pattern[1:], text) < d:
            for x in NUCLEOTIDES:
                neighborhood.add(x + text)
        else:
            neighborhood.add(pattern[0] + text)
    return neighborhood


def neighbors(pattern, d):
    neighborhood = {pattern}
    for j in range(d):
        neighborhood |= immediate_neighbors(pattern_prime)
    return neighborhood


def mutations(pattern, d):
    if d == 0:
        return [pattern]
    elif d == 1:
        return immediate_neighbors(pattern)
    else:
        return set(itertools.chain.from_iterable(
            mutations(mutation, d - 1) for mutation in immediate_neighbors(pattern)
        ))


def neighborhood_brute(pattern, d):
    k = len(pattern)
    return set(
      mutation for mutation in generate_patterns(k)
      if hamming_distance(pattern, mutation) <= d
    )


def neighborhood(pattern, d):
    return mutations(pattern, d)


def neighborhoods(patterns, d):
    return set(itertools.chain.from_iterable(
        neighbors(pattern, d) for pattern in set(patterns)
    ))
