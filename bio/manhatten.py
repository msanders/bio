import numpy as np


def manhatten_tourist(n: int, m: int, down: np.ndarray, right: np.ndarray):
    s = np.zeros([n + 1, m + 1], dtype='i')
    s[0][0] = 0
    for i in range(1, n + 1):
        s[i][0] = s[i - 1][0] + down[i - 1][0]

    for j in range(1, m + 1):
        s[0][j] = s[0][j - 1] + right[0][j - 1]

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s[i][j] = max(s[i - 1][j] + down[i-1][j], s[i][j - 1] + right[i][j-1])
    return s[n][m]


def parse_matrix(text: str) -> np.ndarray:
    return np.asarray([
        [int(x) for x in line.split()] for line in text.strip().splitlines()
    ], dtype='i')


def parse_tourist_input(text: str) -> (int, int, np.ndarray, np.ndarray):
    lines = text.strip().splitlines()
    n, m = lines[0].split()
    n, m = int(n), int(m)
    down = parse_matrix("\n".join(lines[1:][:n]))
    right = parse_matrix("\n".join(lines[2 + n:]))
    return n, m, down, right


def lcs_path(v: str, w: str) -> np.ndarray:
    s = np.zeros((len(v) + 1,len(w) + 1), dtype='i')
    for i in range(len(v)):
        for j in range(len(w)):
            if v[i] == w[j]:
                s[i + 1][j + 1] = s[i][j] + 1
            else:
                s[i + 1][j + 1] = max(s[i + 1][j], s[i][j + 1])
    return s


def lcs_backtrack(s: np.ndarray) -> np.ndarray:
    backtrack = np.zeros([len(s), len(s[0])], dtype='str')

    for i in range(1, len(s)):
        for j in range(1, len(s[0])):
            if s[i][j] == s[i - 1][j]:
                backtrack[i][j] = "↓"
            elif s[i][j] == s[i][j - 1]:
                backtrack[i][j] = "→"
            else:
                backtrack[i][j] = "↘"
    return backtrack


def output_lcs(backtrack: np.ndarray, v: str, w: str):
    lcs = []
    i, j = len(v), len(w)
    while i > 0 and j > 0:
        if backtrack[i][j] == "↓":
            i -= 1
        elif backtrack[i][j] == "→":
            j -= 1
        elif backtrack[i][j] == "↘":
            lcs.append(v[i - 1])
            i -= 1
            j -= 1
    return "".join(reversed(lcs))


def main():
    text = """AACCTTGG
    ACACTGTGA"""
    v, w = text.split()
    s = lcs_path(v, w)
    backtrack = lcs_backtrack(s)
    print(backtrack)
    lcs = output_lcs(backtrack, s, v, w)
    print(lcs)
    return

    text = """
    6 15
    2 0 4 2 0 1 4 4 1 1 2 4 2 3 3 3
    0 3 2 2 0 1 4 1 1 3 1 3 3 4 3 3
    4 1 0 3 1 2 1 0 3 1 2 3 2 0 4 4
    4 1 1 4 2 1 3 4 0 3 4 1 3 3 4 0
    4 1 4 4 0 3 1 0 3 0 3 3 3 0 4 4
    1 0 3 0 0 1 4 4 0 2 2 3 1 2 0 2
    -
    1 3 1 3 4 0 1 2 0 4 1 0 0 2 0
    0 3 1 1 0 2 4 4 4 1 1 2 1 1 3
    3 3 1 1 2 0 4 0 4 1 4 3 2 4 2
    4 0 2 1 0 0 0 1 1 2 4 4 4 3 3
    2 2 3 3 4 2 1 3 2 3 4 3 1 0 1
    0 2 3 3 0 3 4 0 3 4 3 3 2 0 0
    3 4 4 1 0 1 3 0 2 3 2 2 2 0 3
    """
    n, m, down, right = parse_tourist_input(text)
    print(manhatten_tourist(n, m, down, right))


if __name__ == '__main__':
    main()
