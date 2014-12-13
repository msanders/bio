# cython: language_level=3
import numpy as np
cimport numpy as np
import cython

@cython.boundscheck(False)
@cython.wraparound(False)
def kmer_probability(object profile, unicode kmer):
    rows = ['A', 'C', 'G', 'T']
    cdef np.longdouble_t probability = 1.0
    for i, x in enumerate(kmer):
        probability *= profile[rows.index(x)][i]
    return probability
