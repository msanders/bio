from ..composition import string_composition, genome_path_string
from ._lib import find_datasets
import numpy as np


def test_string_composition():
    dataset = [
        (("CAATCCAAC", 5), "AATCC ATCCA CAATC CCAAC TCCAA")
    ] + find_datasets("string_composition", lambda x: x.splitlines())

    for input_sample, output_sample in dataset:
        text, k = input_sample
        k = int(k)
        output = string_composition(text, k)
        output_sample = output_sample.split()
        assert np.array_equal(output, output_sample), "{0} != {1}".format(output, output_sample)


def test_genome_path_string():
    dataset = [
        ("ACCGA CCGAA CGAAG GAAGC AAGCT", "ACCGAAGCT")
    ] + find_datasets("genome_path_string")

    for input_sample, output_sample in dataset:
        kmers = input_sample.split()
        output = genome_path_string(kmers)
        output_sample = output_sample
        assert output == output_sample, "{0} != {1}".format(output, output_sample)


def main():
    test_string_composition()
    test_genome_path_string()
    print("Success!")

if __name__ == '__main__':
    main()