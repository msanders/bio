import numpy as np


def string_composition(text: str, k: int) -> [str]:
    kmercount = len(text) - k + 1
    kmers = np.empty(kmercount, dtype=object)
    for i in range(kmercount):
        kmers[i] = text[i:i + k]
    kmers.sort()
    return kmers


def genome_path_string(kmers: [str]):
    """
    String Spelled by a Genome Path Problem: Reconstruct a string from its
                                             genome path.
    Input: A sequence of k-mers Pattern1, … ,Pattern[n] such that the last k - 1
           symbols of Pattern[i] are equal to the first k-1 symbols of
           Pattern[i+1] for 1 ≤ i ≤ n-1.
    Output: A string Text of length k+n-1 such that the i-th k-mer in Text is
            equal to Pattern[i] (for 1 ≤ i ≤ n).
    """
    k = len(kmers[0])
    pattern = kmers[0]
    for kmer in kmers[1:]:
        pattern += kmer[-1]
    return pattern


def main():
    with open("bio/data/dataset_198_3.txt") as f:
        kmers = f.read().strip().split()
        print(genome_path_string(kmers))
    #print("\n".join(string_composition("CAATCCAAC", 5)))


if __name__ == '__main__':
    main()
