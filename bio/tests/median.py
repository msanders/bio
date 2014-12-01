from ..median import median_string, distance_between_pattern_and_dna
from ._lib import find_datasets

def test_median_string():
    dataset = [
        (("AAATTGACGCAT GACGACCACGTT CGTCAGCGCCTG GCTGAGCACCGG AGTACGGGACAG", 3),
          "GAC"),
    ] + find_datasets("median_string", lambda x: x.splitlines())

    for input_sample, output_sample in dataset:
        dna, k = input_sample
        dna, k = dna.split(), int(k)
        output = median_string(dna, k)
        assert output_sample in output, "{0} != {1}".format(output, output_sample)


def test_distance_between_pattern_and_dna():
    dataset = [
        (("AAA", "TTACCTTAAC GATATCTGTC ACGGCGTTCG CCCTAAAGAG CGTCAGAGGT"), 5),
    ] + find_datasets("distance_between_pattern_and_dna",
                      lambda x: x.splitlines())

    for input_sample, output_sample in dataset:
        pattern, dna = input_sample
        dna = dna.split()
        output = distance_between_pattern_and_dna(pattern, dna)
        output_sample = int(output_sample)
        assert output_sample == output, "{0} != {1}".format(output, output_sample)


def main():
    test_median_string()
    test_distance_between_pattern_and_dna()
    print("Success!")


if __name__ == '__main__':
    main()
