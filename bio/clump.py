from numberpattern import compute_frequencies_dict as compute_frequencies
from numberpattern import number_to_pattern, pattern_to_number
import numberpattern
import sys


def find_clumps_brute(genome, k, l, t):
    clumps = set()
    for i in range(len(genome)):
        if i + l <= len(genome):
            window = genome[i:][:l]
            for key, count in compute_frequencies(window, k).iteritems():
                if count >= t:
                    clumps.add(key)
    return clumps


def baton(i):
    baton = "/-\\|/-\\|"
    return baton[i % len(baton)]


def find_clumps(genome, k, l, t):
    clumps = set()

    frequency_array = numberpattern.compute_frequencies(genome[:l], k)
    for i, value in enumerate(frequency_array):
        if value >= t:
            clumps.add(i)

    for i in range(1, len(genome) - l):
        first_pattern = genome[i - 1:k + i - 1]
        frequency_array[pattern_to_number(first_pattern)] -= 1

        last_pattern = genome[i + l - k: i + l]
        j = pattern_to_number(last_pattern)
        frequency_array[j] += 1

        if frequency_array[j] >= t:
            clumps.add(j)

        if i % 10000 == 0:
            sys.stdout.write("Loading entry {0} {1}\r".format(
                i, baton(i // 10000)
            ))

    return set(number_to_pattern(i, k) for i in clumps)


def main():
    with open("data/E-coli.txt") as f:
        clumps = find_clumps(f.read().strip(), 9, 500, 3)
        print(" ".join(clumps))
        print("{0} unique clumps".format(len(set(clumps))))
    #print(" ".join(find_clumps(
    #    "CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA",
    #    5, 50, 4
    #)))

if __name__ == '__main__':
    main()
    #cProfile.run('main()')
