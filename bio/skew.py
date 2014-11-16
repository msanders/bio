def skew(genome):
    skews = []
    for i in range(len(genome) + 1):
        if i == 0:
            skews.append(0)
        else:
            previous_nucleotide = genome[i - 1]
            if previous_nucleotide == "G":
                skews.append(skews[i - 1] + 1)
            elif previous_nucleotide == "C":
                skews.append(skews[i - 1] - 1)
            else:
                skews.append(skews[i - 1])
        
    return skews


def skew_minimum(genome):
    skews = skew(genome)
    minimum = min(skews)
    return [i for i, value in enumerate(skews) if value == minimum]


def main():
    genome = "TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT"
    print(skew_minimum(genome))
    print(" ".join(map(str, skew(genome))))
    with open("dataset_7_6.txt") as f:
        print(skew_minimum(f.read().strip()))

if __name__ == '__main__':
    main()
