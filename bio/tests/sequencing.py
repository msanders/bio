from ..sequencing import (
    score, cyclopeptide_sequence, cyclopeptide_scored_sequence,
    convolution_cyclopeptide_sequence
)
from ._lib import find_datasets
import numpy as np


def test_cyclic_score():
    dataset = [
        (("NQEL", "0 99 113 114 128 227 257 299 355 356 370 371 484"), 11),
    ] + find_datasets(
        "cyclic_score", lambda x: x.splitlines(), lambda y: int(y)
    )

    for input_sample, output_sample in dataset:
        peptide, spectrum = input_sample
        spectrum = np.asarray([int(x) for x in spectrum.split()], dtype='i')
        output = score(peptide, spectrum, cycle=True)
        assert output == output_sample, "{0} != {1}".format(output, output_sample)


def test_linear_score():
    dataset = [
        (("NQEL", "0 99 113 114 128 227 257 299 355 356 370 371 484"), 8),
    ] + find_datasets(
        "linear_score", lambda x: x.splitlines(), lambda y: int(y)
    )

    for input_sample, output_sample in dataset:
        peptide, spectrum = input_sample
        spectrum = np.asarray([int(x) for x in spectrum.split()], dtype='i')
        output = score(peptide, spectrum, cycle=False)
        assert output == output_sample, "{0} != {1}".format(output, output_sample)


def test_cyclopeptide_sequence():
    dataset = [
        ("0 113 128 186 241 299 314 427", { "186-128-113", "186-113-128",
                                            "128-186-113", "128-113-186",
                                            "113-186-128", "113-128-186" })
    ] + find_datasets("cyclopeptide_sequence",
                      output_transform=lambda x: set(x.split()))

    for input_sample, output_sample in dataset:
        spectrum = np.asarray([int(x) for x in input_sample.split()],
                              dtype='i')
        output = cyclopeptide_sequence(spectrum)
        assert np.array_equal(output, output_sample), "{0} != {1}".format(output, output_sample)


def test_cyclopeptide_scored_sequence():
    dataset = [
        ((10, "0 71 113 129 147 200 218 260 313 331 347 389 460"), {
            "113-147-71-129"
        })
    ]

    for input_sample, output_sample in dataset:
        n, spectrum = input_sample
        spectrum = np.asarray([int(x) for x in spectrum.split()], dtype='i')
        output = cyclopeptide_scored_sequence(n, spectrum)
        assert output_sample.issubset(output), "!{0}.issubset({1})".format(
            output_sample,
            output
        )

def test_convolution_cyclopeptide_sequence():
    dataset = [
        ((20, 60, "57 57 71 99 129 137 170 186 194 208 228 265 285 299 307 323 356 364 394 422 493"), {
            "99-71-137-57-72-57"
        })
    ] + find_datasets("convolution_sequence",
                      input_transform=lambda x: x.rstrip().splitlines(),
                      output_transform=lambda x: set(x.split()))

    for input_sample, output_sample in dataset:
        m, n, spectrum = input_sample
        m, n = int(m), int(n)
        spectrum = np.asarray([int(x) for x in spectrum.split()], dtype='i')
        output = convolution_cyclopeptide_sequence(m, n, spectrum)
        assert output_sample.issubset(output), "!{0}.issubset({1})".format(
            output_sample,
            output
        )


def main():
    test_cyclic_score()
    test_linear_score()
    test_cyclopeptide_sequence()
    test_cyclopeptide_scored_sequence()
    test_convolution_cyclopeptide_sequence()
    print("Success!")

if __name__ == '__main__':
    main()
