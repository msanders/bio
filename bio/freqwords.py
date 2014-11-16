from patcount import patcount


def frequent_words(text, k):
    frequent_patterns = set()
    count = []
    for i in range(len(text) - k):
        pattern = text[i:][:k]
        count.append(patcount(text, pattern))
    maxcount = max(count)
    for i in range(len(text) - k):
        if count[i] == maxcount:
            frequent_patterns.add(text[i:][:k])
    return frequent_patterns


def main():
    text = "ACGTTGCATGTCGCATGATGCATGAGAGCT"
    k = 4
    print(sorted(frequent_words(text, k)))

    with open("dataset_2_9.txt", "r") as f:
        print(frequent_words(f.read(), 12))

if __name__ == '__main__':
    main()
