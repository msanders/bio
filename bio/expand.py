from .mass import MASS_TABLE, linearspectrum, suffix_spectrum
from .csequencing import score, overlapping, trim_scores
import numpy as np

def expand_spectrum(peptide: str,
                    acid: str,
                    old_spectrum: np.ndarray,
                    parent_spectrum: np.ndarray) -> np.ndarray:
    return overlapping(
        parent_spectrum, np.concatenate((
            old_spectrum,
            suffix_spectrum(peptide, acid)
        ))
    )


def expand(peptides: list):
    new_peptides = []
    for acid in MASS_TABLE:
        for peptide in peptides:
            if not peptide:
                new_peptides.append(acid)
            else:
                new_peptides.append(peptide + acid)
    return new_peptides


def expanded(leaderboard: list,
             masses: dict,
             spectrums: dict,
             parent_spectrum: np.ndarray) -> (list, dict, dict):
    expanded_leaderboard = []
    expanded_masses = {}
    expanded_spectrums = {}
    parent_mass = parent_spectrum[-1]

    for acid in MASS_TABLE:
        if acid == 'K' or acid == 'L': # Treat K/Q and L/I as equal.
            continue

        acid_mass = MASS_TABLE[acid]
        for peptide in leaderboard:
            if not peptide:
                if acid_mass <= parent_mass:
                    expanded_leaderboard.append(acid)
                    expanded_masses[acid] = acid_mass
                    expanded_spectrums[acid] =  overlapping(
                        parent_spectrum,
                        np.array((0, acid_mass), dtype='i')
                    )
            else:
                new_mass = masses[peptide] + acid_mass
                if new_mass > parent_mass:
                    continue

                old_spectrum = spectrums[peptide]
                suffix = peptide + acid
                if suffix not in expanded_masses:
                    expanded_leaderboard.append(suffix)
                    expanded_masses[suffix] = new_mass
                    expanded_spectrums[suffix] = expand_spectrum(
                        peptide, acid, old_spectrum, parent_spectrum
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
