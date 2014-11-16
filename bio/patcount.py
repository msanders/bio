def patcount(text, pattern):
    count = 0
    for i in range(len(text) - len(pattern)):
        if text[i:][:len(pattern)] == pattern:
            count += 1
    return count


def main():
    print(patcount(
        "ACAACTATGCATACTATCGGGAACTATCCT",
        "ACTAT"
    ))

    with open("dataset_2_6.txt", "r") as f:
        print(patcount(f.read(), "ATGGTCGAT"))

if __name__ == '__main__':
    main()
