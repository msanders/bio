import math
import itertools
import numpy as np
import numpy.ma as ma
from collections import Counter


def delta(a: [int]) -> [int]:
    delta_a = np.empty(len(a) * len(a), dtype='i')
    i = 0
    for x in a:
        for y in a:
            delta_a[i] = x - y
            i += 1
    delta_a.sort()
    return delta_a


def array_without(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    mask = np.ones(len(a), dtype=bool)
    for x in set(b):
        count = len(np.where(b == x)[0])
        for i in np.where(a == x)[0][:count]:
            mask[i] = False
    return a[mask]


def list_issubset(a: list, b: list) -> bool:
    b1 = list(b)
    for x in a:
        if x not in b1:
            return False
        b1.remove(x)
    return True


def turnpike(delta_a: [int]) -> [int]:
    """
    Turnpike Problem: Given all pairwise distances between points on a line
                      segment, reconstruct the positions of those points.
    Input: A collection of integers L.
    Output: A set of integers A such that ∆A = L.
    => −7 −5 −4 −3 −2 −2 0 0 0 0 2 2 3 4 5 7
       (0-7), (2-7), (0-4), (4-7), (2-4), (0-2), (0-0), (2-2), (4-4), (7-7),
       (2-0), (4-2), (7-4), (4-0), (7-2), (7-0)
    => 0, 2, 4, 7
    """
    delta_a = np.asarray(delta_a, dtype='i')
    a_len = int(math.sqrt(len(delta_a)))
    valid_entries = sorted({x for x in delta_a if x >= 0})

    start = np.asarray([0], dtype='i')
    possible_points = [(start, array_without(delta_a, delta(start)))]
    for i, x in enumerate(valid_entries):
        new_points = []
        for points, delta_points in possible_points:
            if len(points) >= a_len:
                continue

            new_deltas = np.concatenate((points - x, x - points, [0]))
            new_subset = array_without(delta_points, new_deltas)
            if len(new_subset) == len(delta_points) - len(new_deltas):
                appended_points = np.append(points, x)
                new_points.append((
                    appended_points,
                    new_subset
                ))
        possible_points += new_points

    return (x for x, y in possible_points if len(x) == a_len)


def main():
    a = [0, 2, 4, 7, 8]
    print(a)
    print(delta(a))
    for solution in turnpike(delta(a)):
        print(" ".join(map(str, solution)))

    #with open("bio/data/dataset_110_1.txt") as f:
    #    a = [int(x) for x in f.read().strip().split()]
    #    for solution in turnpike(a):
    #        print(" ".join(map(str, solution)))


if __name__ == '__main__':
    main()
