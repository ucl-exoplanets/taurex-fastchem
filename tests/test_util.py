"""Test utilities."""

import numpy as np


def test_abundance_to_string():
    """Test abundance to string conversion."""
    from taurex_fastchem.util import _abundance_to_string

    elements = ["H", "He", "C", "N", "O"]
    abundances = np.array([0.74, 0.24, 0.01, 0.01, 0.01])

    expected = "# Need this in the beginning\nH 7.40000000e-01\nHe 2.40000000e-01\nC 1.00000000e-02\nN 1.00000000e-02\nO 1.00000000e-02\n"

    assert _abundance_to_string(elements, abundances) == expected


def test_get_selected_default():
    """Test whether get_selected_elements returns the default elements."""
    from taurex_fastchem.util import (
        get_selected_default,
        default_elements,
        default_abundances,
    )

    assert get_selected_default(None) == (default_elements, default_abundances)


def test_get_selected_default_choice():
    """Test whether get_selected_elements returns only selected elements."""
    from taurex_fastchem.util import (
        get_selected_default,
        default_elements,
        default_abundances,
    )

    selected_elements = ["H", "He", "O", "C", "N"]

    expected_abundances = [
        default_abundances[default_elements.index(elem)] for elem in selected_elements
    ]

    assert len(selected_elements) == 5

    assert get_selected_default(selected_elements) == (
        selected_elements,
        expected_abundances,
    )


def test_fix_element_abundance_order():
    """Test whether the element abundance order is fixed."""
    from taurex_fastchem.util import fix_element_abundance_order

    elements = ["O", "He", "H", "e-", "C", "N"]
    abundances = [0.01, 0.24, 0.74, 0.05, 0.01, 0.01]

    fixed_elements = ["H", "He", "O", "C", "N", "e-"]
    fixed_abundances = [0.74, 0.24, 0.01, 0.01, 0.01, 0.05]

    assert fix_element_abundance_order(elements, abundances) == (
        fixed_elements,
        fixed_abundances,
    )


def test_fix_element_abundance_order_no_change():
    """Test whether the element abundance order is fixed."""
    from taurex_fastchem.util import fix_element_abundance_order

    elements = ["H", "He", "O", "C", "N"]
    abundances = [0.74, 0.24, 0.01, 0.01, 0.01]

    assert fix_element_abundance_order(elements, abundances) == (
        elements,
        abundances,
    )


def test_create_fastchem_default_args():
    """Test wheter we can create a PyFastchem object"""
    from taurex_fastchem.util import _create_pyfastchem_, default_elements
    from taurex_fastchem import DEFAULT_LOGK_WO_IONS

    fc = _create_pyfastchem_(DEFAULT_LOGK_WO_IONS)

    assert fc is not None

    assert fc.getElementNumber() == len(default_elements)


def test_create_fastcham_filter_element():
    """Test wheter we can create a PyFastchem object, while filtering elements"""
    from taurex_fastchem.util import _create_pyfastchem_, default_elements
    from taurex_fastchem import DEFAULT_LOGK_WO_IONS

    fc = _create_pyfastchem_(DEFAULT_LOGK_WO_IONS, ["H", "He", "O", "C", "N"])

    assert fc is not None

    assert fc.getElementNumber() == 5


def test_fix_fastchem_species():
    """Test whether fastchem species are fixed."""
    from taurex_fastchem.util import fix_fastchem_species

    taurex_species = ["NH3", "H2O", "OH"]
    fastchem_species = ["H3N1", "H2O1", "H1O1", "PH3", "H2O2", "H1O2"]

    new_species = fix_fastchem_species(taurex_species, fastchem_species)

    assert new_species == ["NH3", "H2O", "OH", "PH3", "H2O2", "H1O2"]
