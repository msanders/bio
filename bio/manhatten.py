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


def main():
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
