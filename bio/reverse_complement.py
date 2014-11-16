COMPLEMENTS = {
    "A": "T",
    "G": "C",
    "T": "A",
    "C": "G"
}


def complement(nucleotide):
    return COMPLEMENTS[nucleotide]


def reverse_complement(strand):
    return "".join(reversed([COMPLEMENTS[n] for n in strand]))


def main():
    with open("dataset_3_2.txt") as f:
        print(reverse_complement(f.read().strip()))

if __name__ == '__main__':
    main()
