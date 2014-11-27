from ..median import median_string
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


def main():
    test_median_string()
    print("Success!")


if __name__ == '__main__':
    main()
