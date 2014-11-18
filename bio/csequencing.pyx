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


def score(object peptide, object spectrum, cycle=None):
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
def trim_against(object ar1, object ar2, int n):
    cdef:
        np.intp_t i, j
        np.intp_t[:] sorted_indexes
        list trimmed

    if len(ar2) == 0:
        return ar1
    elif n == 0:
        return []

    sorted_indexes = np.argsort(ar2)
    trimmed = []

    for i from len(ar2) > i >= 0:
        j = sorted_indexes[i]
        if i <= len(ar2) - n:
            break

        trimmed.append(ar1[j])

    # Add equivalent items in original order.
    if len(ar2) >= n:
        last_value = ar2[sorted_indexes[<np.intp_t>len(ar2) - n]]
        for i from 0 <= i <= len(ar2) - n:
            j = sorted_indexes[i]
            v = ar2[j]
            if v != last_value:
                continue

            trimmed.append(ar1[j])

    return trimmed
