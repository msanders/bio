from .mutations import neighborhood
from .hamming import hamcount, hamming_distance
from functools import reduce
import operator
import numpy as np

def motif_enumeration(dna: [str], k: int, d: int) -> {str}:
    """
    Input: A collection of strings Dna, and integers k and d.
    Output: All (k, d)-motifs in Dna.
    """
    neighborhoods = []
    patterns = set()
    for text in dna:
        for start in (text[i:] for i in range(k)):
            for i in range(len(start) - k):
                pattern = start[i:i + k]
                for pattern_prime in neighborhood(pattern, d + 1):
                    found = True
                    for pat in dna:
                        if hamcount(pat, pattern_prime, d) == 0:
                            found = False
                            break
                    if found:
                        patterns.add(pattern_prime)
        return patterns


def kmer_probability(profile: [float], kmer: str) -> float:
    rows = ['A', 'C', 'G', 'T']
    return reduce(operator.mul, (profile[rows.index(x)][i]
                  for i, x in enumerate(kmer)))

def profile_most_probable(profile: [float], text: str, k: int) -> str:
    """
    Input: A string Text, an integer k, and a 4 × k matrix Profile.
    Output: A Profile-most probable k-mer in Text.
    """
    max_probability = -1
    most_probable = None
    for start in (text[i:] for i in range(k)):
        for i in range(len(start) - k):
            kmer = start[i:i + k]
            probability = kmer_probability(profile, kmer)
            if probability > max_probability:
                most_probable = kmer
                max_probability = probability
    return most_probable


def profile_from_dna(dna: [str]) -> [float]:
    profile = np.zeros((4, len(dna[0])))
    rows = ['A', 'C', 'G', 'T']
    for text in dna:
        for i, nucleotide in enumerate(text):
            row = rows.index(nucleotide)
            profile[row][i] += 1.0
    profile /= float(len(dna))
    return profile


def parse_profile(input: str) -> [float]:
    return np.array([[np.float(x) for x in row.split()]
                     for row in input.strip().splitlines()])


def find_consensus(motifs: [str]) -> str:
    profile = profile_from_dna(motifs)
    rows = ['A', 'C', 'G', 'T']
    consensus = [""] * len(motifs[0])
    max_probability = [0] * len(motifs[0])
    for i, row in enumerate(profile):
        for j, column in enumerate(row):
            if column > max_probability[j]:
                consensus[j] = rows[i]
                max_probability[j] = column
    return "".join(consensus)


def score(motifs: [str]) -> int:
    consensus = find_consensus(motifs)
    score = 0
    for motif in motifs:
        score += hamming_distance(consensus, motif)
    return score


def greedy_motif_search(dna: [str], k: int, t: int):
    best_motifs = [text[:k] for text in dna]
    best_score = score(best_motifs)
    #print(best_motifs)
    #print(best_score)
    for i in range(len(dna[0]) - k):
        kmer = dna[0][i:i + k]
        motif = [kmer]
        for j in range(1, min(t, len(dna))):
            #print("Motif: {0}".format(motif))
            profile = profile_from_dna(motif)
            most_probable = profile_most_probable(profile, dna[j], k)
            if most_probable:
                motif.append(most_probable)

        motif_score = score(motif)
        if motif_score < best_score:
            best_motifs = motif
            best_score = motif_score

    return best_motifs


def main():
    profile = parse_profile(
        """
        0.258 0.197 0.318 0.242 0.273 0.273 0.258 0.318 0.273 0.212 0.409 0.258 0.212 0.212 0.258
        0.227 0.258 0.242 0.197 0.242 0.167 0.303 0.152 0.242 0.242 0.182 0.303 0.273 0.288 0.288
        0.242 0.227 0.227 0.273 0.227 0.212 0.152 0.273 0.197 0.273 0.197 0.242 0.258 0.273 0.242
        0.273 0.318 0.212 0.288 0.258 0.348 0.288 0.258 0.288 0.273 0.212 0.197 0.258 0.227 0.212
        """
    )
    #print(profile)
    print(profile_most_probable(profile, "AGTGACCATTTACCCGAGGCTTCCGACAAGTAGGATTTACTCTTGGGATCCTGAGCTCGTAGTGTTTGAGGCCCTCTAGGTGCTTTTGAGGGATTGCCTGACGTGTCGGGATAAAAGCAATCAGCCACTGGACGTCTCTTGCTACCCATTTATGGCCAAATGAGTGACCAGCTATAAAGTGGTGGTGCCTTTACGAATAAAGAGTGTTCGAGCGAAAGCAAATTCGACCTGGAAGGGAGTTAGGAATTTGACGAATTGCAGACCAGTTAGAGCCACCGTGCTTCTGCTTGCCGCGGAATGGAATCACAATACAGGAGGGTCACATTTTCCCTCCAAAGACGACAGTTGCAGAGTATCCCGCGATCTACAGCGTTTGCTGTGGCCACGATCTGAATAGGGTCAGTGGTTAATTCTCCATACTATCGATGAAATCTGTGTAAATCAACAGAGGGGCGGCGACTGCCACGTCATAGAATCAGTTAAACTAGGTGCGATCCGGCCTGTAGTGCGACAATATAGCTTACTATCTTCATATGGCAACCAGCAAATTAGATTTCGGTCACAGTAACCCGATCTGTTACGCGAGGTGCACTCAGATGAGTCGCTCGGGAAACCCATAGCAACGAGGATTAGCGAATTGAAGTACCTAGGGAGGAACCCGCGAACGATAGTTTTAAGTCTGCGCCGAATCGTTTTACTGAACGTGTATCGGTGCCAATATGCTGAATGAACACACTACACGTTTATAGCCCAAAACGTTTAAAGCGAGGTGTATCTTGGAACCGGCAAGCATACGCTGCACATCTAGCTTCGCATATAACGGCAGTTGAGGAATCTAATGGGACACTTGCGGTGTGCTTTGTCGGAGAAGTTAGAAACAGCAGTTGGTTCGGCGGGGTCCTGACTTAGGGCACGGCGAACATACAGCGCTAGGAGTCAGCAAGCATATTGGAAACTTAACGCCTTCTTGAACCCTTAATCTAAATGGCAAAAATATGCG", 15))
    #print(" ".join(motif_enumeration("GCCTTGACGGTCAGGCTTACTTGAA CGTGAATGACGTCTTGTTCTGGATG CCGACCAGGTGTCTTAGAGAGATGT CTATGCCGGCGATGTGCCTTCCTGC GTCTTCAGCTCGGGTGTGGAGGGGC CAAGGTGCTGCAAAACGTTGGTCTT".split(), 5, 1)))

if __name__ == '__main__':
    main()
