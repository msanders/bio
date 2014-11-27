import math

def entropy(*args):
    return -sum(0 if x == 0 else x * math.log(x, 2) for x in args)

def main():
    entropies = [
        entropy(0.2, 0.1, 0.0, 0.7),
        entropy(0.2, 0.6, 0.2, 0.0),
        entropy(0.0, 0.0, 1.0, 0.0),
        entropy(0.0, 0.0, 1.0, 0.0),
        entropy(0.0, 0.0, 0.9, 0.1),
        entropy(0.0, 0.0, 0.9, 0.1),
        entropy(0.9, 0.0, 0.1, 0.0),
        entropy(0.1, 0.4, 0.0, 0.5),
        entropy(0.1, 0.1, 0.0, 0.8),
        entropy(0.1, 0.2, 0.0, 0.7),
        entropy(0.3, 0.4, 0.0, 0.3),
        entropy(0.0, 0.6, 0.0, 0.4),
    ]

    print(entropies)
    print(sum(entropies))

if __name__ == '__main__':
    main()
