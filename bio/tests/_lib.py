import os
import random
from ..mass import MASS_TABLE

DATA_DIR = os.path.join("bio", "data" ,"tests")
AMINOS = [x for x in MASS_TABLE.keys() if x != "K" and x != "L"]


def read(path: str) -> str:
    with open(path) as f:
        return f.read()


def find_datasets(name, input_transform=None, output_transform=None):
    if input_transform is None:
        input_transform = lambda x: x
    if output_transform is None:
        output_transform = lambda x: x

    datasets = []
    for fname in os.listdir(DATA_DIR):
        input_path = os.path.join(DATA_DIR, fname)
        output_path = os.path.join(DATA_DIR, fname.replace("_input", "_output"))
        if fname.startswith(name + "_input") and os.path.exists(output_path):
            datasets.append((
                input_transform(read(input_path).strip()),
                output_transform(read(output_path).strip())
            ))
    return datasets


def random_peptide(length: int) -> str:
    return "".join(random_amino() for _ in range(length))


def random_amino() -> str:
    return random.choice(AMINOS)


def spectrum_difference(ar1: [int], ar2: [int]) -> [int]:
    difference = list(ar1)
    for x in ar2:
        if x in difference:
            difference.remove(x)
    return difference
