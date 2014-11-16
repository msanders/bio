from protein_translation import codon_table, chunks
from reverse_complement import reverse_complement


def all_chunks(string, chunksize):
    result = []
    for i in range(chunksize):
        result += list(chunks(string[i:], chunksize))
    return result


def encode_peptide(text, peptide):
    """
    Input: A DNA string Text, an amino acid string Peptide, and the array
           GeneticCode.
    Output: All substrings of Text encoding Peptide (if any such substrings
            exist).
    """
    chunklen = 3
    substrings = all_chunks(text, chunklen * len(peptide))

    output = []
    for text in substrings:
        for pattern in [text, reverse_complement(text)]:
            match = True
            for i, c in enumerate(peptide):
                if codon_table(pattern[chunklen * i:][:chunklen]) != c:
                    match = False
                    break
                else:
                    a, b = pattern[:chunklen], pattern[chunklen:][:chunklen]
                    #print("{0}, {1} => {2}, {3} => {4}".format(
                    #    pattern, a, codon_table(a), b, codon_table(b)
                    #))
            if match:
                output.append(text)
    return output


def main():
    #print(" ".join(encode_peptide(
    #    "ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA",
    #    "MA"
    #)))

    #with open("data/peptide_encoding_in.txt") as f:
    #    encoded = encode_peptide(f.read().strip(), "KEVFEPHYY")
    #    #print("ENCODED")
    #    print("\n".join(sorted(set(encoded))))

    #    #with open("data/peptide_encoding_out.txt") as f2:
    #    #    print("EXPECTED")
    #    #    expected = f2.read().strip().split()
    #    #    print("\n".join(sorted(set(expected))))
    #    #    print(set(expected) == set(encoded))
    #    #    print("False negatives: {0}".format(set(expected) - set(encoded)))
    #    #    print("False positives: {0}".format(set(encoded) - set(expected)))

    #with open("data/dataset_96_8.txt") as f:
    #    encoded = encode_peptide(f.read().strip(), "LVWWRERE")
    #    print("\n".join(encoded))

    with open("data/B_Brevis.txt") as f:
        data = "".join(line.rstrip() for line in f)
        encoded = encode_peptide(data, "VKLFPWFNQY")
        print("\n".join(encoded))

if __name__ == '__main__':
    main()
