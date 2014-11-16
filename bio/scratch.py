from .mass import MASS_TABLE, linearspectrum, cyclospectrum, mass, subpeptides
import numpy as np

def putl(lst):
    print(list(lst))


def spectrum_append1(peptide, amino, spectrum):
    start = np.array((0, MASS_TABLE[peptide[0]]))
    return np.sort(
        np.concatenate((
            start,
            start + MASS_TABLE[amino],
            start + MASS_TABLE[peptide[-1]],
            [spectrum[-1] + MASS_TABLE[amino]]
        ))
    )


def spectrum_append(peptide, amino, spectrum):
    start = np.array([0] + [MASS_TABLE[x] for x in peptide[:-1]])
    return np.sort(
        np.concatenate((
            start,
            start + MASS_TABLE[amino],
            #start + MASS_TABLE[peptide[-1]],
            [spectrum[-1] + MASS_TABLE[amino]]
        ))
    )


def added_peptides(peptide, amino):
    old_subpeptides = subpeptides(peptide, False)
    new_subpeptides = subpeptides(peptide + amino, False)
    changed_peptides = set(new_subpeptides) - set(old_subpeptides)
    added_masses = set(linearspectrum(peptide + amino)) - set(linearspectrum(peptide))
    return {x: mass(x) for x in changed_peptides if mass(x) in added_masses}


def new_masses_naive(peptide, amino, spectrum):
    amino_mass = MASS_TABLE[amino]
    return amino_mass + np.fromiter(
        (mass(peptide[:k]) for k in range(len(peptide) + 1)),
        dtype='i'
    )

def new_masses(peptide, amino, spectrum):
    amino_mass = MASS_TABLE[amino]
    masses = np.empty(len(peptide) + 1, dtype='i')
    masses[0] = amino_mass
    for i in range(1, len(peptide) + 1):
        j = len(peptide) - i
        masses[i] = masses[i - 1] + MASS_TABLE[peptide[j]]
    return masses


def test_prepend(peptide, amino):
    old_spectrum = linearspectrum(peptide)
    new_peptide = peptide + amino
    new_spectrum = linearspectrum(new_peptide)
    print("{0} => {1} (peptide mass: {2})".format(peptide, new_peptide, MASS_TABLE[amino]))
    print("Previous spectrum: {0}".format(old_spectrum))
    print("New spectrum: {0}".format(new_spectrum))
    print("Added elements: {0}".format(sorted(set(new_spectrum) - set(old_spectrum))))
    #print("Added peptides: {0}".format(added_peptides(peptide, amino)))
    #print("Cyclic peptides: {0}".format(set(cyclospectrum(new_peptide)) - set(linearspectrum(new_peptide))))
    print("Estimated: {0}".format(sorted(set(new_masses(peptide, amino, old_spectrum)))))


af = "AF"
daf = "DAF"
cdaf = "CDAF"
ycdaf = "YCDAF"
mycdaf = "MYCDAF"
lmycdaf = "LMYCDAF"
a_mass = MASS_TABLE["A"]
d_mass = MASS_TABLE["D"]
f_mass = MASS_TABLE["F"]
c_mass = MASS_TABLE["C"]
y_mass = MASS_TABLE["Y"]

print("a, d, f, c, y = {0}".format((a_mass, d_mass, f_mass, c_mass, y_mass)))
test_prepend(af, "D"); print()
test_prepend(daf, "C"); print()
test_prepend(cdaf, "Y"); print()
test_prepend(ycdaf, "M"); print()
test_prepend(mycdaf, "L")
