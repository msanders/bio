import numpy as np


def string_composition(text: str, k: int) -> [str]:
    kmercount = len(text) - k + 1
    kmers = np.empty(kmercount, dtype=object)
    for i in range(kmercount):
        kmers[i] = text[i:i + k]
    kmers.sort()
    return kmers


def genome_path_string(kmers: [str]) -> str:
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


def paired_composition(text: str, k: int, d: int) -> [(str, str)]:
    kmercount = len(text) - k + 1
    pairs = []
    for i in range(kmercount):
        d_offset = i + d + k
        if d_offset + k <= len(text):
            pairs.append((text[i:i + k], text[d_offset:d_offset + k]))
    pairs.sort()
    return pairs


def string_spelled_by_gapped_patterns(patterns: [(str, str)], d: int) -> str:
    first_patterns = [x[0] for x in patterns]
    second_patterns = [x[1] for x in patterns]
    k = len(first_patterns[0])
    prefix_string = genome_path_string(first_patterns)
    suffix_string = genome_path_string(second_patterns)
    for i in range(k + d + 1, len(prefix_string)):
        if prefix_string[i] != suffix_string[i - k - d]:
            return None
    return prefix_string + suffix_string[-(k + d):]


def parse_patterns(text: str) -> [(str, str)]:
    text = text.replace("(", "").replace(")", "")
    return [x.split("|") for x in text.split()]


def main():
    with open("bio/data/dataset_6206_7.txt") as f:
        patterns = parse_patterns(f.read().strip())
        #patterns = parse_patterns("GACC|GCGC ACCG|CGCC CCGA|GCCG CGAG|CCGG GAGC|CGGA")
        print(string_spelled_by_gapped_patterns(patterns, 200))
    #print(paired_composition("TAATGCCATGGGATGTT", 3, 2))
    #print(" ".join("(" + "|".join(x) + ")" for x in paired_composition("TAATGCCATGGGATGTT", 3, 2)))
    #with open("bio/data/dataset_198_3.txt") as f:
    #    kmers = f.read().strip().split()
    #    print(genome_path_string(kmers))
    #print("\n".join(string_composition("CAATCCAAC", 5)))


if __name__ == '__main__':
    main()
