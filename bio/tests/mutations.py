from ..mutations import neighborhood
from ._lib import find_datasets


def test_neighborhood():
    dataset = [
        (("ACG", 1), "CCG TCG GCG AAG ATG AGG ACA ACC ACT ACG")
    ] + find_datasets("neighborhood", lambda x: x.splitlines())

    for input_sample, output_sample in dataset:
        pattern, d = input_sample
        d = int(d)
        output_sample = set(output_sample.split())
        output = neighborhood(pattern, d)
        assert output == output_sample, "{0} != {1}".format(output, output_sample)


def main():
    test_neighborhood()
    print("Success!")

if __name__ == '__main__':
    main()

