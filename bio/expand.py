from .mass import (
    AMINO_MASSES, EXTENDED_ALPHABET, MASS_TABLE,
    linearspectrum, suffix_spectrum
)
from .csequencing import score, overlapping, trim_scores
import numpy as np

def expand_spectrum(peptide: tuple,
                    acid_mass: int,
                    old_spectrum: np.ndarray,
                    parent_spectrum: np.ndarray) -> np.ndarray:
    return overlapping(
        parent_spectrum, np.concatenate((
            old_spectrum,
            suffix_spectrum(peptide, acid_mass)
        ))
    )


def expand(peptides: list):
    new_peptides = []
    for acid in MASS_TABLE:
        if acid == 'K' or acid == 'L': # Treat K/Q and L/I as equal.
            continue

        for peptide in peptides:
            if not peptide:
                new_peptides.append(acid)
            else:
                new_peptides.append(peptide + acid)
    return new_peptides


def expanded(leaderboard: list,
             masses: dict,
             spectrums: dict,
             parent_spectrum: np.ndarray,
             use_extended_alphabet=None) -> (list, dict, dict):
    if use_extended_alphabet is None:
        use_extended_alphabet = False
    expanded_leaderboard = []
    expanded_masses = {}
    expanded_spectrums = {}
    parent_mass = parent_spectrum[-1]
    alphabet = EXTENDED_ALPHABET if use_extended_alphabet else AMINO_MASSES

    for acid_mass in alphabet:
        for peptide in leaderboard:
            if not peptide:
                if acid_mass <= parent_mass:
                    leader = (acid_mass,)
                    expanded_leaderboard.append(leader)
                    expanded_masses[leader] = acid_mass
                    expanded_spectrums[leader] =  overlapping(
                        parent_spectrum,
                        np.array((0, acid_mass), dtype='i')
                    )
            else:
                new_mass = masses[peptide] + acid_mass
                if new_mass > parent_mass:
                    continue

                old_spectrum = spectrums[peptide]
                suffix = peptide + (acid_mass,)
                if suffix not in expanded_masses:
                    expanded_leaderboard.append(suffix)
                    expanded_masses[suffix] = new_mass
                    expanded_spectrums[suffix] = expand_spectrum(
                        peptide, acid_mass, old_spectrum, parent_spectrum
                    )

    return (expanded_leaderboard, expanded_masses, expanded_spectrums)


def trim(leaderboard: list, parent_spectrum: np.ndarray, n: int) -> list:
    scores = np.fromiter(
        (score(x, parent_spectrum, False) for x in leaderboard),
        dtype='i',
        count=len(leaderboard)
    )

    return trim_scores(leaderboard, scores, n)


def trim_parallel(leaderboard: list, spectrum: np.ndarray, n: int) -> list:
    with ProcessPoolExecutor() as executor:
        scores = executor.map(ScoreCall(spectrum, False), leaderboard)

    return trim_scores(
        leaderboard,
        np.fromiter(scores, dtype='i', count=len(leaderboard)),
        n
    )


def trim_expanded_spectrums(leaderboard, spectrums, comparison_spectrum, n):
    scores = np.fromiter(
        (len(spectrums[x]) for x in leaderboard),
        dtype='i',
        count=len(leaderboard)
    )

    return trim_scores(leaderboard, scores, n)


# Workaround for Python's annoying lack of closure support for concurrency.
class ScoreCall(object):
    def __init__(self, spectrum: np.ndarray, cycle: bool):
        self.spectrum = spectrum
        self.cycle = cycle

    def __call__(self, peptide: str):
        return score(peptide, self.spectrum, self.cycle)
