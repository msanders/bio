import functools
import sys


@functools.lru_cache(2 ** 10)
def min_num_coins(money: int, coins: [int]) -> int:
    if money == 0:
        return 0

    min_coins = float("inf")
    for coin in coins:
        if money >= coin:
            num_coins = min_num_coins(money - coin, coins)
            if num_coins + 1 < min_coins:
                min_coins = num_coins + 1
    return min_coins


class Memoize(object):
    def __init__(self, f):
        self.f = f
        self.memo = {}

    def __call__(self, *args, **kwargs):
        if (args, str(kwargs)) in self.memo:
            return self.memo[(args, str(kwargs))]

        out = self.f(*args, **kwargs)
        self.memo[(args, str(kwargs))] = out
        return out


def main():
    #memo = Memoize(min_num_coins)
    coins = (17,5,3,1)
    sys.setrecursionlimit(15000)
    print(min_num_coins(19526, coins))
    #values = [memo(i, coins, memo) for i in range(13, 23)]
    #print(" ".join([ str(x) for x in values ]))
    #for i in range(13, 23):
    #    print("MinNumCoins({0}): {1}".format(i, memo(i, coins, memo)))
    #print(Memoize(min_num_coins)(50, (5, 4, 1)))


if __name__ == '__main__':
    main()
