def patmatch(text, pattern):
    indexes = []
    for i in range(len(text) - len(pattern)):
        if text[i:][:len(pattern)] == pattern:
            indexes.append(i)
    return indexes


def main():
    with open("Vibrio_cholerae.txt") as f:
        print(" ".join(map(str, patmatch(
            f.read().strip(),
            "CTTGATCAT"
        ))))

if __name__ == '__main__':
    main()

