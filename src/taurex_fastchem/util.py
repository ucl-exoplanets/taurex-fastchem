"""Utility functions for FastChem."""

import typing as t
import numpy as np
import numpy.typing as npt
import pathlib
import pyfastchem


def _abundance_to_string(
    elements: t.List[str], abundances: npt.NDArray[np.float64]
) -> str:
    """Build abundance string for FastChem.

    Args:
        elements (List[str]): List of element symbols.
        abundances (np.float64): List of element abundances.

    Returns:
        str: Abundance string for FastChem.

    """
    # Build the abundance string
    abundance_string = "# Need this in the beginning\n"

    if len(elements) != len(abundances):
        raise ValueError("Length of elements and abundances must be equal.")

    for elem, abund in zip(elements, abundances):
        abundance_string += f"{elem} {abund:.8e}\n"

    return abundance_string


# Default elements and abundances
default_elements = [
    "H",
    "He",
    "O",
    "Al",
    "Ar",
    "C",
    "Ca",
    "Cl",
    "Co",
    "Cr",
    "Cu",
    "F",
    "Fe",
    "Ge",
    "K",
    "Mg",
    "Mn",
    "N",
    "Na",
    "Ne",
    "Ni",
    "P",
    "S",
    "Si",
    "Ti",
    "V",
    "Zn",
    "e-",
]
default_abundances = [
    12.00,
    10.93,
    8.69,
    6.45,
    6.40,
    8.43,
    6.34,
    5.50,
    4.99,
    5.64,
    4.19,
    4.56,
    7.50,
    3.65,
    5.03,
    7.60,
    5.43,
    7.83,
    6.24,
    7.93,
    6.22,
    5.41,
    7.12,
    7.51,
    4.95,
    3.93,
    4.56,
    13.1139,
]


def fix_element_abundance_order(
    elements: t.List[str], abundances: t.List[float]
) -> t.Tuple[t.List[str], t.List[float]]:
    """Fix the element abundance order to ensure H, He and O are
    in that order."""
    # Ensure H, He and O are first three elements
    if "O" in elements:
        idx = elements.index("O")
        elements.insert(0, elements.pop(idx))
        abundances.insert(0, abundances.pop(idx))

    if "He" in elements:
        idx = elements.index("He")
        elements.insert(0, elements.pop(idx))
        abundances.insert(0, abundances.pop(idx))

    if "H" in elements:
        idx = elements.index("H")
        elements.insert(0, elements.pop(idx))
        abundances.insert(0, abundances.pop(idx))

    # Ensure e- is last element
    if "e-" in elements:
        idx = elements.index("e-")
        elements.append(elements.pop(idx))
        abundances.append(abundances.pop(idx))

    return elements, abundances


def get_selected_default(
    elements: t.Optional[str] = None,
    element_choice: t.Optional[str] = None,
    abundance_choice: t.Optional[str] = None,
) -> t.Tuple[t.List[str], t.List[float]]:
    """Get selected default elements and abundances."""
    element_choice = element_choice or default_elements
    abundance_choice = abundance_choice or default_abundances

    if elements is None:
        return default_elements, default_abundances

    selected_elements = []
    selected_abundances = []

    for elem in elements:
        if elem in default_elements:
            idx = default_elements.index(elem)
            selected_elements.append(elem)
            selected_abundances.append(default_abundances[idx])

    selected_elements, selected_abundances = fix_element_abundance_order(
        selected_elements, selected_abundances
    )

    return selected_elements, selected_abundances


def load_abundance_file(filename: t.Union[str, pathlib.Path]):
    """Load the abundance file."""
    filename = pathlib.Path(filename)

    if not filename.exists():
        raise FileNotFoundError(f"File {filename} not found")

    with open(filename, "r") as f:
        lines = f.readlines()

    elements = []
    abundances = []

    for line in lines:
        if line.startswith("#"):
            continue
        elem, abund = line.split()
        elements.append(elem)
        abundances.append(float(abund))

    # Ensure H, He and O are first three elements

    elements, abundances = fix_element_abundance_order(elements, abundances)

    return elements, abundances


def _create_pyfastchem_(
    log_k: pathlib.Path,
    selected_elements: t.Optional[t.List[str]] = None,
    element_file: t.Optional[pathlib.Path] = None,
) -> str:
    """Create the PyFastchem object.

    Here we will generate a PyFastchem object with the selected elements.
    We wil then go on to modify the elements and their abundances later on
    in the code. We jsut need fastchem to load it in.

    """
    import tempfile
    import os

    # Load the elements
    elements, abundances = default_elements, default_abundances
    if element_file is not None:
        elements, abundances = load_abundance_file(element_file)

    # Get the selected elements
    elements, abundances = get_selected_default(selected_elements, elements, abundances)

    # Create the abundance string
    abundance_string = _abundance_to_string(elements, abundances)
    # Now create
    element_file = tempfile.NamedTemporaryFile(mode="w", delete=False)
    element_file.write(abundance_string)
    elem_filename = element_file.name
    element_file.close()

    # Create the fastchem object
    fc = pyfastchem.FastChem(str(elem_filename), str(log_k), 1)

    # Remove the temporary file
    os.unlink(elem_filename)

    return fc


def element_count(mol):
    """Count the number of elements in a molecule."""
    import re

    s = re.findall("([A-Z][a-z]?)([0-9]*)", mol)
    compoundweight = 0
    elems_dict = {}
    for element, count in s:
        if count == "":
            count = 1
        if element not in elems_dict:
            elems_dict[element] = int(count)
        else:
            elems_dict[element] += count
    return elems_dict


def fix_fastchem_species(
    taurex_species: t.List[str], fastchem_species: t.List[str]
) -> t.List[str]:
    """Fix Fastchem species to match TauREx names."""
    new_species_list = []
    for s in fastchem_species:
        gas_name = s
        elem_check = element_count(s)

        for av in taurex_species:
            elem_t_check = element_count(av)
            if elem_check == elem_t_check:
                gas_name = av
                break
        new_species_list.append(gas_name)
    return new_species_list
