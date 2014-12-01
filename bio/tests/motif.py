from ..motif import (
    motif_enumeration, parse_profile, profile_most_probable,
    greedy_motif_search, sampling_iterator, randomized_motif_search,
    gibbs_sampler
)
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
        assert output_sample == output, "{0} != {1}".format(output, output_sample)


def test_greedy_motif_search():
    def parse_input(text):
        lines = text.splitlines()
        return "\n".join(lines[:-2]), lines[-2], lines[-1]

    dataset = [
        (("GGCGTTCAGGCA AAGAATCAGTCA CAAGGAGTTCGC CACGTCAATCAC CAATAATATTCG", 3, 5),
          "CAG CAG CAA CAA CAA")
    ] + find_datasets("greedy_motif_search", parse_input)

    for input_sample, output_sample in dataset:
        dna, k, t = input_sample
        dna, k, t = dna.split(), int(k), int(t)
        output = greedy_motif_search(dna, k, t)
        output_sample = output_sample.split()
        assert output_sample == output, "{0} != {1}".format(output, output_sample)


def test_greedy_pseudocount_motif_search():
    def parse_input(text):
        lines = text.splitlines()
        return "\n".join(lines[:-2]), lines[-2], lines[-1]

    dataset = [
        (("GGCGTTCAGGCA AAGAATCAGTCA CAAGGAGTTCGC CACGTCAATCAC CAATAATATTCG", 3, 5),
          "TTC ATC TTC ATC TTC")
    ] + find_datasets("greedy_psuedocount_motif_search", parse_input)

    for input_sample, output_sample in dataset:
        dna, k, t = input_sample
        dna, k, t = dna.split(), int(k), int(t)
        output = greedy_motif_search(dna, k, t, pseudocount=True)
        output_sample = output_sample.split()
        assert output_sample == output, "{0} != {1}".format(output, output_sample)


def test_randomized_motif_search():
    def parse_input(text):
        lines = text.splitlines()
        return "\n".join(lines[:-2]), lines[-2], lines[-1]

    dataset = [
        (("CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG "
          "TAGTACCGAGACCGAAAGAAGTATACAGGCGT TAGATCAAGTTTCAGGTGCACGTCGGTGAACC "
          "AATCCACCAGCTCCACGTGCAATGTTGGCCTA", 8, 5),
          "TCTCGGGG CCAAGGTG TACAGGCG TTCAGGTG TCCACGTG")
    ] + find_datasets("randomized_motif_search", parse_input)

    for input_sample, output_sample in dataset:
        dna, k, t = input_sample
        dna, k, t = dna.split(), int(k), int(t)
        output = sampling_iterator(randomized_motif_search, dna, k, t,
                                   concurrent=True, pseudocount=True)
        output_sample = output_sample.split()
        assert output_sample == output, "{0} != {1}".format(output, output_sample)
        output_sample = output_sample.split()
        assert output_sample == output, "{0} != {1}".format(output, output_sample)


def main():
    test_motif_enumeration()
    test_profile_most_probable()
    test_greedy_motif_search()
    test_greedy_pseudocount_motif_search()
    test_randomized_motif_search()
    print("Success!")

if __name__ == '__main__':
    main()
