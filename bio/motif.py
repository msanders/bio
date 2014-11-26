from .mutations import neighborhood
from .hamming import hamcount

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


def main():
    print(" ".join(motif_enumeration("GCCTTGACGGTCAGGCTTACTTGAA CGTGAATGACGTCTTGTTCTGGATG CCGACCAGGTGTCTTAGAGAGATGT CTATGCCGGCGATGTGCCTTCCTGC GTCTTCAGCTCGGGTGTGGAGGGGC CAAGGTGCTGCAAAACGTTGGTCTT".split(), 5, 1))) 

if __name__ == '__main__':
    main()
