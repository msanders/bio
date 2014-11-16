# cython: language_level=3
from .cmass import cyclospectrum, linearspectrum
from .cmass cimport _MASS_TABLE as MASS_TABLE
from cython import nogil, gil
from cython.parallel import prange
from concurrent.futures import ProcessPoolExecutor
import cython
cimport numpy as np
import numpy as np


@cython.boundscheck(False)
@cython.wraparound(False)
cdef np.ndarray[np.uint8_t, ndim=1] overlapping_indices(np.ndarray ar1, 
                                                        np.ndarray ar2):
    cdef:
        np.ndarray[np.uint8_t] found = np.zeros(len(ar1), dtype=np.uint8)
        np.intp_t i
        int x, y

    for x in ar2:
        for i, y in enumerate(ar1):
            if x == y and found[i] == 0:
                found[i] = 1
                break

    return found


@cython.boundscheck(False)
@cython.wraparound(False)
def overlapping(np.ndarray ar1, np.ndarray ar2):
    return ar1[overlapping_indices(ar1, ar2).nonzero()]


def score(str peptide, object spectrum, cycle=None):
    """
    Input: An amino acid string Peptide and a collection of integers Spectrum.
    Output: The score of Peptide against Spectrum, Score(Peptide, Spectrum).
    """
    cycle = cycle if cycle is not None else True
    return np.count_nonzero(overlapping_indices(
        spectrum,
        cyclospectrum(peptide) if cycle else linearspectrum(peptide)
    ))


@cython.boundscheck(False)
@cython.wraparound(False)
def trim_scores(list leaderboard, np.ndarray scores, int n):
    cdef:
        np.intp_t i, j
        np.intp_t[:] sorted_indexes
        int last_score = 0
        int max_score
        list trimmed

    if scores.size == 0:
        return leaderboard

    sorted_indexes = np.argsort(scores)
    trimmed = []

    # Add scores in descending order.
    for i from scores.size > i >= 0:
        j = sorted_indexes[i]
        v = scores[j]
        if i <= scores.size - n:
            break

        trimmed.append(leaderboard[j])

    # Add equivalent items in original order.
    if scores.size >= n:
        last_score = scores[sorted_indexes[<int>scores.size - n]]
        for i from 0 <= i <= scores.size - n:
            j = sorted_indexes[i]
            v = scores[j]
            if v != last_score:
                continue

            trimmed.append(leaderboard[j])

    #max_score = scores[sorted_indexes[<int>scores.size - 1]]
    #print("Score: {0}, {1}".format(max_score, last_score))
    return trimmed
