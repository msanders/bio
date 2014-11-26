from ..motif import motif_enumeration
from ._lib import find_datasets


def test_motif_enumeration():
    dataset = [
        (("ATTTGGC TGCCTTA CGGTATC GAAAATT", 3, 1), "ATA ATT GTT TTT")
    ] + find_datasets("motif_enumeration", lambda x: x.splitlines())

    for input_sample, output_sample in dataset:
        dna, k, d = input_sample
        dna, k, d = dna.split(), int(k), int(d)
        output_sample = set(output_sample.strip().split())
        output = motif_enumeration(dna, k, d)
        assert output == output_sample, "{0} != {1}".format(output, output_sample)


def main():
    test_motif_enumeration()
    print("Success!")

if __name__ == '__main__':
    main()
