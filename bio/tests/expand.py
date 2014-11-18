from ..expand import trim, trim_against
from ._lib import find_datasets, random_peptide
import numpy as np
import random


def test_trim_against():
    for length in range(1, 20):
        n = random.randint(0, length)
        peptide = random_peptide(length)
        order = np.random.rand(length)
        sorted_order = list(reversed(sorted(order)))

        trimmed_peptide = trim_against(peptide, order, n)
        assert len(trimmed_peptide) >= n
        assert len(trimmed_peptide) == n or (
            trimmed_peptide[-1] == trimmed_peptide[-2]
        )

        order = list(order)
        for i in range(len(trimmed_peptide)):
            j = order.index(sorted_order[i])
            assert trimmed_peptide[i] == peptide[j]


def test_trim():
    dataset = [
        (("LAST ALST TLLT TQAS", "0 71 87 101 113 158 184 188 259 271 372", 2),
          "LAST ALST"),
    ] + find_datasets("trim", lambda x: x.splitlines())

    for input_sample, output_sample in dataset:
        leaderboard, spectrum, n = input_sample
        leaderboard = leaderboard.split()
        spectrum = np.asarray([int(x) for x in spectrum.split()], dtype='i')
        n = int(n)
        output_sample = output_sample.split()
        output = trim(leaderboard, spectrum, n)
        assert np.array_equal(output, output_sample), "{0} != {1}".format(output, output_sample)


def main():
    test_trim_against()
    test_trim()
    print("Success!")

if __name__ == '__main__':
    main()
