from hamming import approx_frequent_words, approx_frequent_complements, approx_frequent_complements_brute_parallel
from skew import skew_minimum
import itertools

def approx_dnaa(genome, k, d, l):
    skews = skew_minimum(genome)
    position = skews[0]
    print("Found skews: {0}".format(skews))
    print("Looking for most frequent {0}-mers (with {1} mismatch and "
          "reverse complements) within a window of {2} starting "
          "at {3}".format(k, d, l, position))

    window = genome[position:position + l]
    complements = approx_frequent_complements(window, k, d)
    found_complements = {x for x in complements if x in window}
    missing_complements = complements - found_complements
    print("Found {0} complements ({1} not in genome): {2}".format(
        len(complements), len(missing_complements), missing_complements
    ))

    return found_complements

def read_fasta(f):
    return "".join(line.rstrip() for line in itertools.islice(f, 1, None))

def main():
    #with open("data/E-coli.txt", "r") as f:
    #    contents = f.read().strip()
    #    print(approx_dnaa(contents, 9, 1, 500))

    with open("data/Salmonella_enterica.txt", "r") as f:
        contents = read_fasta(f)
        #print(approx_dnaa(contents, 9, 0, 500))
        #print(approx_dnaa(contents, 9, 1, 500))
        print(approx_dnaa(contents, 9, 1, 1000))
        print(approx_dnaa(contents, 9, 1, 2000))


if __name__ == '__main__':
    main()
