# -*- coding: utf-8 -*-
import pyximport; pyximport.install()  # pylint: disable=F0401
from chamming import hamming_distance, hamcount  # pylint: disable=F0401
from concurrent.futures import ProcessPoolExecutor
from mutations import neighborhoods
from numberpattern import all_nucleotides
from reverse_complement import reverse_complement
import sys


# Workaround for Python's annoying lack of closure support for concurrency.
class HamcountCalc(object):
    def __init__(self, text, d):
        self.text = text
        self.d = d

    def __call__(self, pattern):
        return hamcount(self.text, pattern, self.d)


def approx_patmatch(text, pattern, d):
    indexes = []
    for i in range(len(text) - len(pattern) + 1):
        if hamming_distance(text[i:i + len(pattern)], pattern) <= d:
            indexes.append(i)
    return indexes


def approx_frequent_words_brute(text, k, d):
    patterns = list(all_nucleotides(k, stringify=True))
    counts = []

    print("{0} patterns".format(len(patterns)))
    for i, pattern in enumerate(patterns):
        if i % 10000 == 0:
            sys.stdout.write("Counting pattern {0}\r".format(i))
        counts.append(hamcount(text, pattern, d))

    maxcount = max(counts)
    return set(
        pattern for pattern, count in zip(patterns, counts)
        if count == maxcount
    )

def approx_frequent_words_brute_parallel(text, k, d):
    patterns = list(all_nucleotides(k, stringify=True))
    counts = []

    print("{0} patterns".format(len(patterns)))
    with ProcessPoolExecutor() as executor:
        counts = executor.map(HamcountCalc(text, d), patterns)

    counts = list(counts)
    maxcount = max(counts)
    return set(
        pattern for pattern, count in zip(patterns, counts)
        if count == maxcount
    )


def approx_frequent_complements_brute_parallel(text, k, d):
    patterns = list(all_nucleotides(k, stringify=True))
    print("{0} patterns".format(len(patterns)))
    counts = []
    with ProcessPoolExecutor() as executor:
        counts = executor.map(HamcountCalc(text, d), patterns)

    patterns = dict(zip(patterns, counts))
    maxcount = 0
    frequent_complements = set()
    for pattern, count in patterns.items():
        if pattern not in frequent_complements:
            complement = reverse_complement(pattern)
            sumcount = count + patterns[complement]
            if sumcount > maxcount:
                maxcount = sumcount
                frequent_complements = {pattern, complement}
            elif sumcount == maxcount:
                frequent_complements |= {pattern, complement}

    return frequent_complements


def approx_frequent_words(text, k, d):
    """
    Input: A string Text as well as integers k and d.
           (You may assume k ≤ 12 and d ≤ 3.)
    Output: All most frequent k-mers with up to d mismatches in Text.
    """
    patterns = neighborhoods((text[i:i + k] for i in range(len(text) - k )), d)
    print("{0} patterns".format(len(patterns)))

    counts = [hamcount(text, pattern, d) for pattern in patterns]
    maxcount = max(counts)
    return set(
        pattern for pattern, count in zip(patterns, counts)
        if count == maxcount
    )

def approx_frequent_complements(text, k, d):
    """
    Input: A DNA string Text as well as integers k and d.
    Output: All k-mers Pattern maximizing the sum Count_d(Text, Pattern) +
            Count_d(Text, PatternComplement) over all possible k-mers.
    """
    patterns = neighborhoods((text[i:i + k] for i in range(len(text) - k )), d)
    patterns |= set(reverse_complement(pattern) for pattern in patterns)

    print("{0} patterns".format(len(patterns)))
    counts = [hamcount(text, pattern, d) for pattern in patterns]

    counts = dict(zip(patterns, counts))
    maxcount = 0
    frequent_complements = set()
    for pattern, count in counts.items():
        if pattern not in frequent_complements:
            complement = reverse_complement(pattern)
            sumcount = count + counts[complement]
            if sumcount > maxcount:
                maxcount = sumcount
                frequent_complements = {pattern, complement}
            elif sumcount == maxcount:
                frequent_complements |= {pattern, complement}

    print("maxcount: {0}".format(maxcount))
    return frequent_complements

def main():
    #print(approx_frequent_words_brute_parallel("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4, 1))
    #with open("data/dataset_9_8.txt")  as f:
    #    contents = f.read().strip()
    #    print(" ".join(sorted(approx_frequent_words(contents, 8, 3))))
    #print(" ".join(sorted(approx_frequent_complements("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4, 1))))
    with open("data/frequent_words_mismatch_data.txt") as f:
        contents = f.read().strip()
        print(" ".join(sorted(approx_frequent_complements(contents, 10, 2))))
    #with open("data/dataset_9_7.txt") as f:
    #    contents = f.read().strip()
    #    #print(hamcount(contents , "GCACACAGAC", 2))
    #    #print(hamcount(contents, "GCGCACACAC", 2))
    #    print(" ".join(sorted(approx_frequent_words(contents, 9, 3))))
    #with open("dataset_9_7.txt") as f:
    #    contents = f.read().strip()
    #    #print(hamcount(contents, "GAAATTAAAT", 2))
    #    #print(hamcount(contents, "TTGATGATGA", 2))
    #    #print(hamcount(contents, "AAATGGTAAA", 2))
    #    #print(" ".join(sorted(approx_frequent_words(contents, 10, 2))))
    #    print(" ".join(sorted(approx_frequent_words_brute_parallel(contents, 10, 2))))

if __name__ == '__main__':
    main()
