CODON_TABLE = None

def codon_table(peptide, fallback=None):
    global CODON_TABLE
    if CODON_TABLE is None:
        CODON_TABLE = {}
        with open("data/RNA_codon_table_1.txt") as f:
            for line in f:
                pair = line.split()
                if len(pair) > 1:
                    codon, amino = pair
                    CODON_TABLE[codon] = amino
    return CODON_TABLE.get(peptide.replace("T", "U"), fallback)


def chunks(string, chunksize):
    return (string[i:chunksize + i] for i in range(0, len(string), chunksize) 
            if chunksize + i <= len(string))


def translate_protein(pattern):
    """
    Input: An RNA string Pattern and the array GeneticCode.
    Output: The translation of Pattern into an amino acid string Peptide.
    """
    return "".join(codon_table(x, "") for x in chunks(pattern, 3))


def main():
    print(translate_protein("AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA"))
    with open("data/dataset_96_5.txt") as f:
        print(translate_protein(f.read()))

if __name__ == '__main__':
    main()
