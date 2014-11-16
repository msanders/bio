# -*- coding: utf-8 -*-
import numpy as np
import pyximport
pyximport.install(setup_args={
    "include_dirs": np.get_include(),
})  # pylint: disable=F0401

from concurrent.futures import ProcessPoolExecutor
from .mass import MASS_TABLE, mass, cyclospectrum, linearspectrum
from .csequencing import score
from .expand import expand, expanded, trim_expanded_spectrums, trim
from itertools import permutations


def peptide_code(peptide: object):
    if isinstance(peptide, str):
        return "-".join(str(MASS_TABLE[c]) for c in peptide)
    code = "-".join(str(x) for x in peptide)
    return code


def cyclopeptide_sequence(spectrum: [int]) -> {str}:
    spectrum = np.asarray(spectrum, dtype='i')
    peptides = [""]
    parent_mass = spectrum[-1]
    found = []
    while peptides:
        peptides = expand(peptides)
        for peptide in peptides.copy():
            if mass(peptide) == parent_mass:
                if np.array_equal(cyclospectrum(peptide), spectrum):
                    found.append(peptide)
                peptides.remove(peptide)
            elif not set(linearspectrum(peptide)).issubset(set(spectrum)):
                peptides.remove(peptide)

    return {peptide_code(x) for x in found}


def cyclopeptide_scored_sequence_naive(n: int, spectrum: [int]) -> {str}:
    """
    Input: A collection of integers Spectrum.
    Output: A cyclic peptide Peptide maximizing Score(Peptide, Spectrum) over
            all peptides Peptide with mass equal to ParentMass(Spectrum).
    """
    spectrum = np.array(spectrum, dtype='i')
    leaderboard = [""]
    leader_peptides = set()
    leader_score = score("", spectrum)
    parent_mass = spectrum[-1]
    while leaderboard:
        leaderboard = expand(leaderboard)
        dead_peptides = []
        print("len: {0}".format(len(leaderboard[0])))
        print("Before: {0}".format(len(leaderboard)))
        for peptide in leaderboard:
            peptide_mass = mass(peptide)
            if peptide_mass == parent_mass:
                peptide_score = score(peptide, spectrum)
                if peptide_score > leader_score:
                    leaderboard = [x for x in leaderboard
                                   if x not in leader_peptides]
                    leader_peptides = {peptide}
                    leader_score = peptide_score
                elif peptide_score == leader_score:
                    leader_peptides.add(peptide)
            elif peptide_mass > parent_mass:
                dead_peptides.append(peptide)

        for peptide in dead_peptides:
            leaderboard.remove(peptide)
        leaderboard = trim(leaderboard, spectrum, n)
        print("After: {0}".format(len(leaderboard)))

    return {peptide_code(x) for x in leader_peptides}


def cyclopeptide_scored_sequence(n: int, spectrum: [int]) -> {str}:
    """
    Input: A collection of integers Spectrum.
    Output: A cyclic peptide Peptide maximizing Score(Peptide, Spectrum) over
            all peptides Peptide with mass equal to ParentMass(Spectrum).
    """
    spectrum = np.asarray(spectrum, dtype='i')
    leaderboard = [""]
    masses = {"": 0}
    spectrums = {"": np.empty(0, dtype='i')}
    leader_peptides = [""]
    leader_score = score("", spectrum)
    parent_mass = spectrum[-1]

    while len(leaderboard) > 0:
        start_count = len(leaderboard)
        leaderboard, masses, spectrums = expanded(
            leaderboard, masses, spectrums, spectrum
        )
        expand_count = len(leaderboard)
        for peptide in leaderboard:
            if masses[peptide] == parent_mass:
                peptide_score = score(peptide, spectrum)
                if peptide_score > leader_score:
                    leader_peptides = [peptide]
                    leader_score = peptide_score
                elif peptide_score == leader_score:
                    leader_peptides.append(peptide)

        leaderboard = trim_expanded_spectrums(leaderboard, spectrums, spectrum, n)
        print("start expand-remove trim: {0} {1} {2}".format(
            start_count, expand_count, len(leaderboard)
        ))

    print("Got {0} leaders with score {1}".format(
        len(leader_peptides),
        leader_score)
    )

    return [peptide_code(x) for x in leader_peptides]


def main():
    #with open("bio/data/Tyrocidine_B1_Spectrum_10.txt") as f:
    #    codes = cyclopeptide_scored_sequence(
    #        1000,
    #        [int(x) for x in f.read().strip().split()]
    #    )

    #    expected = "97-129-97-147-99-71-186-71-113-163-115-71-113-128-103-87-128-101-137-163-114"
    #    print(codes) # Should be 148 peptides of length 11
    #    print(expected)

    #with open("bio/data/leaderboard.txt") as f:
    #    codes = cyclopeptide_scored_sequence(
    #        325,
    #        [int(x) for x in f.read().strip().split()]
    #    )

    #    expected = "97-129-97-147-99-71-186-71-113-163-115-71-113-128-103-87-128-101-137-163-114"
    #    print(codes)
    #    print(expected)

    #with open("bio/data/dataset_102_7.txt") as f:
    #    codes = cyclopeptide_scored_sequence(
    #        319,
    #        [int(x) for x in f.read().strip().split()]
    #    )

    #    print(codes)

    with open("bio/data/dataset_103_1.txt") as f:
        codes = cyclopeptide_scored_sequence(
            1000,
            [int(x) for x in f.read().strip().split()]
        )

        print("\n".join(codes))
        lengths = [len(x[0]) for x in codes]
        for length in sorted(set(lengths)):
            print("{0} of length {1}".format(lengths.count(length), length))

        with open("/Users/msanders/out.txt", "w") as f:
            f.write(" ".join(x[1] for x in codes))

if __name__ == '__main__':
    main()
