from numberpattern import generate_kmers

def prob(alphabet, pattern, t, k):
    count = 0
    matches = 0
    for item in generate_kmers(alphabet, k):
        if pattern in item:
            matches += 1
        count += 1
    print("{0} / {1} matches".format(matches, count))
    return matches / count


def main():
    print(list(generate_kmers("01", 4)))
    print(prob("01", "01", 1, 25))

if __name__ == '__main__':
    main()
