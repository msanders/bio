from ..overlapping import overlapping_patterns
from ._lib import find_datasets


def test_overlapping_patterns():
    def parse_sample_output(text: str) -> dict:
        output = text.replace(" -> ", ":").split()
        out_dict = {}
        for x in output:
            key, value = x.split(":")
            if key not in out_dict:
                out_dict[key] = []
            out_dict[key].append(value)
        return out_dict

    dataset = [
        ("ATGCG GCATG CATGC AGGCA GGCAT", 
         "AGGCA -> GGCAT CATGC -> ATGCG GCATG -> CATGC GGCAT -> GCATG")
    ] + find_datasets("overlapping_patterns")

    for input_sample, output_sample in dataset:
        patterns = input_sample.split()
        output_sample = parse_sample_output(output_sample)
        output = overlapping_patterns(patterns)
        assert output == output_sample, "{0} != {1}".format(output, output_sample)


def main():
    test_overlapping_patterns()
    print("Success!")

if __name__ == '__main__':
    main()
