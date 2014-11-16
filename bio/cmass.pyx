# cython: language_level=3
# distutils: extra_compile_args = -Wno-unused-function
import numpy as np
cimport numpy as np
import cython

ctypedef int[:] intarray 
cdef np.ndarray _EXPANDED_ALPHABET = np.arange(57, 201, dtype='i')
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
EXPANDED_ALPHABET = _EXPANDED_ALPHABET


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


@cython.boundscheck(False)
@cython.wraparound(False)
cdef np.ndarray spectrum(str peptide, bint cyclic):
    cdef:
        np.intp_t i, j
        int[:] prefix_mass = np.empty(len(peptide) + 1, dtype='i')
    prefix_mass[0] = 0
    for i in range(1, len(peptide) + 1):
        amino = peptide[i - 1]
        prefix_mass[i] = prefix_mass[i - 1] + _MASS_TABLE[amino]

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


def cyclospectrum(str peptide):
    return spectrum(peptide, True)


def linearspectrum(str peptide):
    return spectrum(peptide, False)


def suffix_spectrum(str peptide, str acid):
    cdef:
        np.ndarray masses
        np.intp_t i, j
    masses = np.empty(len(peptide) + 1, dtype='i')
    masses[0] = _MASS_TABLE[acid]
    for i in range(1, len(peptide) + 1):
        j = len(peptide) - i
        masses[i] = masses[i - 1] + _MASS_TABLE[peptide[j]]
    return masses
