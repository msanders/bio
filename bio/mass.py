# -*- coding: utf-8 -*-
import numpy as np
import pyximport
import functools
pyximport.install(setup_args={
    "include_dirs": np.get_include(),
})  # pylint: disable=F0401
from .cmass import (
    mass, cyclospectrum, linearspectrum, suffix_spectrum, convolution,
    MASS_TABLE, EXTENDED_ALPHABET, AMINO_MASSES
)  # pylint: disable=F0401


def subpeptides(peptide, cycle=None):
    if cycle is None:
        cycle = True
    subpeptides = []
    for length in range(1, len(peptide)):
        for i in range(len(peptide)):
            subpeptide = peptide[i:][:length]
            if len(subpeptide) < length:
                if not cycle:
                    continue
                subpeptide += peptide[:length - len(subpeptide)]
            subpeptides.append(subpeptide)
    return subpeptides


def cyclospectrum_brute(peptide, cycle=None):
    return sorted(mass(x) for x in [""] + subpeptides(peptide, cycle) + [peptide])


def linearspectrum_brute(peptide):
    return cyclospectrum_brute(peptide, cycle=False)


def peptide_string(peptide: tuple) -> str:
    alphabet = list(MASS_TABLE.keys())
    values = list(MASS_TABLE.values())
    output = ""
    for num in peptide:
        if num in values:
            output += alphabet[values.index(num)]
        else:
            output += str(num)
    return output


MASSES = {57, 71, 87, 97, 99, 101, 103, 113, 114, 115, 128, 129, 131, 137, 147, 156, 163, 186}
@functools.lru_cache(2 ** 10)
def peptide_count_with_mass(mass: int) -> int:
    count = 1 if mass in MASSES else 0
    count += sum(peptide_count_with_mass(mass - x) for x in MASSES if x < mass)
    return count


def main():
    print(peptide_count_with_mass(1024))
    #print(" ".join(subpeptides("ELEL")))
    #print(" ".join([""] + subpeptides("NQEL") + ["NQEL"]))
    print(" ".join(map(str, linearspectrum_brute("NQEL"))))
    print(" ".join(map(str, list(linearspectrum("NQEL")))))

    print(" ".join(map(str, cyclospectrum_brute("LEQN"))))
    print(" ".join(map(str, list(cyclospectrum("LEQN")))))
    #print(" ".join(map(str, linear_spectrum("NQEL"))))
    #with open("data/LinearSpectrum.txt") as f:
    #    print(spectrum)
    #    print(spectrum == f.read().strip())

if __name__ == '__main__':
    main()
