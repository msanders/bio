import math
from .hamming import hamming_distance
from .numberpattern import all_nucleotides, number_to_pattern

def min_distance(pattern, text):
    """
    Given a k-mer Pattern and a longer string Text, we use d(Pattern, Text) to
    denote the minimum Hamming distance between Pattern and any k-mer in Text,
    """
    k = len(pattern)
    result = float("inf")
    for i in range(len(text) - k):
        splice = text[i:i + k]
        distance = hamming_distance(pattern, splice)
        if distance < result:
            result = distance
    return result


def min_dna_distance(pattern, dna):
    """
    Given a k-mer Pattern and a set of strings Dna = {Dna1, … , Dnat}, we define
    d(Pattern, Dna) as the sum of distances between Pattern and all strings 
    in Dna.
    """
    return sum(min_distance(pattern, text) for text in dna)


def distance_between_pattern_and_dna(pattern: str, dna: [str]) -> int:
    """
    Input: A string Pattern followed by a collection of strings Dna.
    Output: d(Pattern, Dna).
    """
    k = len(pattern)
    distance = 0
    for text in dna:
        min_distance = float("inf")
        for i in range(len(text) - k + 1):
            kmer = text[i:i + k]
            kmer_distance = hamming_distance(pattern, kmer)
            if min_distance > kmer_distance:
                min_distance = kmer_distance
        distance += min_distance
    return distance


def median_string(dna, k):
    distance = float("inf")
    medians = []
    for i in range(4 ** k):
        pattern = number_to_pattern(i, k)
        pattern_distance = distance_between_pattern_and_dna(pattern, dna)
        if distance > pattern_distance:
            distance = pattern_distance
            medians = [pattern]
        elif distance == pattern_distance:
            medians.append(pattern)
    return medians


def main():
    assert(min_distance("GATTCTCA", "GCAAAGACGCTGACCAA") == 3)
    assert(min_dna_distance("AAA", "TTACCTTAAC GATATCTGTC ACGGCGTTCG CCCTAAAGAG CGTCAGAGGT".split()) == 5)
    print(median_string("AAATTGACGCAT GACGACCACGTT CGTCAGCGCCTG GCTGAGCACCGG AGTACGGGACAG".split(), 3))

if __name__ == '__main__':
    main()
