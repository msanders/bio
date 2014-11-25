# cython: language_level=3
# distutils: extra_compile_args = -Wno-unused-function
import numpy as np
cimport numpy as np
import cython

ctypedef int[:] intarray 
cdef np.ndarray _EXTENDED_ALPHABET = np.arange(57, 201, dtype='i')
cdef np.ndarray _AMINO_MASSES
cdef dict _MASS_TABLE = {
    'V': 99,
    'F': 147,
    'N': 114,
    'E': 129,
    'S': 87,
    'H': 137,
    'I': 113,
    'C': 103,
    'A': 71,
    'P': 97,
    'R': 156,
    'D': 115,
    'L': 113,
    'W': 186,
    'K': 128,
    'G': 57,
    'T': 101,
    'M': 131,
    'Y': 163,
    'Q': 128
}

MASS_TABLE = _MASS_TABLE
AMINO_MASSES = np.unique(np.fromiter(_MASS_TABLE.values(), dtype='i'))
EXTENDED_ALPHABET = _EXTENDED_ALPHABET


def mass(str peptide):
    return sum(_MASS_TABLE[amino] for amino in peptide)


cdef inline int spectrum_size(int n, bint cyclic):
    cdef int size
    if cyclic:
        return n * (n - 1) + 2
    elif n < 2:
        return 2
    else:
        size = n + 2
        while n > 2:
            n -= 1
            size += n
        return size

cdef inline tuple peptide_table(str peptide):
    return tuple(_MASS_TABLE[x] for x in peptide)


@cython.boundscheck(False)
@cython.wraparound(False)
cdef np.ndarray spectrum(tuple peptide, bint cyclic):
    cdef:
        np.intp_t i, j
        int[:] prefix_mass = np.empty(len(peptide) + 1, dtype='i')
    prefix_mass[0] = 0
    for i in range(1, len(peptide) + 1):
        amino_mass = peptide[i - 1]
        prefix_mass[i] = prefix_mass[i - 1] + amino_mass

    cdef:
        int peptide_mass = prefix_mass[len(peptide)]
        np.intp_t k = 1
        np.ndarray output = np.empty(
            spectrum_size(len(peptide), cyclic), dtype='i'
        )

    output[0] = 0
    for i in range(len(peptide)):
        for j in range(i + 1, len(peptide) + 1):
            output[k] = np.int32(prefix_mass[j] - prefix_mass[i])
            k += 1
            if cyclic and i > 0 and j < len(peptide):
                output[k] = np.int32(
                    peptide_mass - (prefix_mass[j] - prefix_mass[i])
                )
                k += 1

    output.sort()
    return output


def cyclospectrum(object peptide):
    peptide = peptide_table(peptide) if isinstance(peptide, str) else peptide
    return spectrum(peptide, True)


def linearspectrum(object peptide):
    peptide = peptide_table(peptide) if isinstance(peptide, str) else peptide
    return spectrum(peptide, False)


def convolution(np.ndarray spectrum):
    """
    Input: A collection of integers Spectrum.
    Output: The list of elements in the convolution of Spectrum. If an element
            has multiplicity k, it should appear exactly k times; you may 
            return the elements in any order.
    """
    differences = []
    for x in spectrum:
        for y in spectrum:
            if x - y > 0:
                differences.append(x - y)
    return differences


def suffix_spectrum(tuple peptide, int acid_mass):
    cdef:
        np.ndarray masses
        np.intp_t i, j
    masses = np.empty(len(peptide) + 1, dtype='i')
    masses[0] = acid_mass
    for i in range(1, len(peptide) + 1):
        j = len(peptide) - i
        masses[i] = masses[i - 1] + peptide[j]
    return masses
