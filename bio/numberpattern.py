from itertools import product
import array

NUCLEOTIDES = "ACGT"

def generate_kmers(alphabet, k, stringify=True):
    output = product(alphabet, repeat=k)
    return output if not stringify else ("".join(x) for x in output)


def all_nucleotides(k, stringify=None):
    return generate_kmers(NUCLEOTIDES, k, stringify)


def pattern_to_number_brute(pattern):
    return list("".join(x) for x in all_nucleotides(
        len(pattern), stringify=True)
    ).index(pattern)


def number_to_pattern_brute(index, k):
    return "".join(list(all_nucleotides(k, stringify=True))[index])


def pattern_to_number_recursive(pattern):
    if not pattern:
        return 0
    nucleotides = "ACGT"
    return 4 * pattern_to_number(pattern[:-1]) + nucleotides.index(pattern[-1])


def number_to_pattern_recursive(index, k):
    if k == 1:
        return NUCLEOTIDES[index]
    quotient, remainder = divmod(index, 4)
    return number_to_pattern(quotient, k - 1) + NUCLEOTIDES[remainder]


def pattern_to_number(pattern):
    num = 0
    for c in pattern:
        num *= 4
        num += NUCLEOTIDES.index(c)
    return num


def number_to_pattern(index, k):
    pattern = []
    quotient = index
    while k > 0:
        quotient, remainder = divmod(quotient, 4)
        pattern.append(NUCLEOTIDES[remainder])
        k -= 1

    return "".join(reversed(pattern))


def compute_frequencies(text, k):
    print("{0} entries".format(4 ** k))
    frequency_array = array.array('h', [0]) * (4 ** k)
    for i in range(len(text) - k):
        pattern = text[i:][:k]
        j = pattern_to_number(pattern)
        frequency_array[j] += 1
    return frequency_array


def compute_frequencies_dict(text, k):
    frequencies = {}
    for i in range(len(text) - k):
        pattern = text[i:][:k]
        frequencies[pattern] = frequencies.get(pattern, 0) + 1
    return frequencies


def faster_frequent_words(text, k):
    frequent_patterns = set()
    frequency_array = compute_frequencies(text, k)
    maxcount = max(frequency_array)
    for i in xrange((4 ** k) - 1):
        if frequency_array[i] == maxcount:
            pattern = number_to_pattern(i, k)
            frequent_patterns.add(pattern)
    return frequent_patterns


def frequency_by_sorting(text, k):
    frequent_patterns = set()
    indexes = []
    count = []
    for i in range(len(text) - k):
        pattern = text[:i][k:]
        indexes.append(pattern_to_number(pattern))
        count.append(1)

    sorted_indexes = sorted(indexes)
    for i in range(1, len(text) - k):
        if sorted_indexes[i] == sorted_indexes[i - 1]:
            count[i] = count[i - 1] + 1

    maxcount = max(count)
    for i in range(len(text) - k):
        if count[i] == maxcount:
            pattern = number_to_pattern(sorted_indexes[i], k)
            frequent_patterns.add(pattern)
    return frequent_patterns


def main():
    print(pattern_to_number("ATGCAA"))
    print(number_to_pattern(5437, 8))
    print(number_to_pattern(2, 1))
    print(frequency_by_sorting("AAGCAAAGGTGGG", 2))
    #with open("dataset_2994_5.txt", "r") as f:
    #    print(" ".join(map(str, compute_frequencies(f.read(), 6))))

if __name__ == '__main__':
    main()
