from ..mass import linearspectrum, cyclospectrum, suffix_spectrum
from ._lib import (
    find_datasets, random_peptide, random_amino, spectrum_difference
)
import numpy as np

def test_linearspectrum():
    dataset = [
        ("NQEL", "0 113 114 128 129 242 242 257 370 371 484")
    ] + find_datasets("linearspectrum")

    for input_sample, output_sample in dataset:
        output_sample = [int(x) for x in output_sample.split()]
        output = linearspectrum(input_sample)
        assert np.array_equal(output, output_sample), "{0} != {1}".format(output, output_sample)


def test_cyclospectrum():
    dataset = [
        ("LEQN", "0 113 114 128 129 227 242 242 257 355 356 370 371 484")
    ] + find_datasets("cyclospectrum")

    for input_sample, output_sample in dataset:
        output_sample = [int(x) for x in output_sample.split()]
        output = cyclospectrum(input_sample)
        assert np.array_equal(output, output_sample), "{0} != {1}".format(output, output_sample)


def test_suffix_spectrum():
    for length in range(1, 20):
        peptide, amino = random_peptide(length), random_amino()
        new_peptide = peptide + (amino,)

        old_spectrum = linearspectrum(peptide)
        new_spectrum = linearspectrum(new_peptide)

        difference = spectrum_difference(new_spectrum, old_spectrum)
        output = suffix_spectrum(peptide, amino)
        assert np.array_equal(difference, output), "{0} != {1}".format(difference, output)


def main():
    test_linearspectrum()
    test_cyclospectrum()
    test_suffix_spectrum()
    print("Success!")

if __name__ == '__main__':
    main()
