#cython: language_level=3, boundscheck=False
def hamming_distance(unicode p, unicode q):
    cdef size_t i
    cdef size_t distance
    distance = 0
    for i in range(len(p)):
        distance += p[i] != q[i]
    return distance


def hamcount(unicode text, unicode pattern, size_t d):
    cdef size_t i
    cdef size_t count
    count = 0
    for i in range(len(text) - len(pattern) + 1):
        if hamming_distance(text[i:i + len(pattern)], pattern) <= d:
            count += 1
    return count
