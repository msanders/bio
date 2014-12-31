from ..manhatten import manhatten_tourist, parse_tourist_input
from ._lib import find_datasets

def test_manhatten_tourist():
    dataset = [(
        """
        4 4
        1 0 2 4 3
        4 6 5 2 1
        4 4 5 2 1
        5 6 8 5 3
        -
        3 2 4 0
        3 2 4 2
        0 7 3 3
        3 3 0 2
        1 3 2 2
        """
    , 34)] + find_datasets("manhatten_tourist")

    for input_sample, output_sample in dataset:
        n, m, down, right = parse_tourist_input(input_sample)
        assert manhatten_tourist(n, m, down, right) == int(output_sample)


def main():
    test_manhatten_tourist()
    print("Success!")


if __name__ == '__main__':
    main()
