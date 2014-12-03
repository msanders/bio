import numpy as np


def string_composition(text: str, k: int) -> [str]:
    kmercount = len(text) - k + 1
    kmers = np.empty(kmercount, dtype=object)
    for i in range(kmercount):
        kmers[i] = text[i:i + k]
    kmers.sort()
    return kmers


def main():
    print("\n".join(string_composition("CAATCCAAC", 5)))


if __name__ == '__main__':
    main()
