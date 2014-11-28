from ..motif import motif_enumeration, parse_profile, profile_most_probable
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


def test_profile_most_probable():
    def parse_input(text):
        lines = text.splitlines()
        return "\n".join(lines[:-2]), lines[-2], lines[-1]

    dataset = [
        (("""
          0.2 0.2 0.3 0.2 0.3
          0.4 0.3 0.1 0.5 0.1
          0.3 0.3 0.5 0.2 0.4
          0.1 0.2 0.1 0.1 0.2
          """, "ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT", 5),
          "CCGAG")
    ] + find_datasets("profile_most_probable", parse_input)

    for input_sample, output_sample in dataset:
        profile, text, k = input_sample
        profile, k = parse_profile(profile), int(k)
        output = profile_most_probable(profile, text, k)
        assert output_sample in output, "{0} != {1}".format(output, output_sample)


def main():
    test_motif_enumeration()
    test_profile_most_probable()
    print("Success!")

if __name__ == '__main__':
    main()
