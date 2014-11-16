# -*- coding: utf-8 -*-
import numpy as np
import pyximport
pyximport.install(setup_args={ 
    "include_dirs": np.get_include(),
})  # pylint: disable=F0401
from .cmass import (
    mass, cyclospectrum, linearspectrum, suffix_spectrum, MASS_TABLE
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


def main():
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
