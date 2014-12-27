from ..alignment import min_num_coins
from ._lib import find_datasets
import numpy as np
import sys

def test_min_num_coins():
    def parse_input(text):
        lines = text.splitlines()
        return int(lines[0]), tuple(int(x) for x in lines[1].split(","))

    dataset = [
        ((40, (50,25,20,10,5,1)), 2)
    ] + find_datasets("min_num_coins", lambda x: parse_input(x))

    for input_sample, output_sample in dataset:
        money, coins = input_sample
        output_sample = int(output_sample)
        output = min_num_coins(money, coins)
        assert output == output_sample, "{0} != {1}".format(output, output_sample)


def main():
    test_min_num_coins()


if __name__ == '__main__':
    sys.setrecursionlimit(15000)
    main()
    print("Success!")
