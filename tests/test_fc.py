"""Test Fastchem class."""

from taurex_fastchem import FastChem
import pytest


def test_fastchem_init():
    """Test fastchem init."""
    fc = FastChem()

    assert len(fc.gases) > 0


def test_fastchem_selected_elements():
    """Test fastchem init."""
    fc_all = FastChem()
    fc = FastChem(selected_elements=["H", "He", "O", "C", "N"])

    assert len(fc.gases) > 0

    assert len(fc_all.gases) > len(fc.gases)


def test_fastchem_various_logk():
    """Test Fastchem with different log k"""
    fc_default = FastChem()
    fc_w_ions = FastChem(with_ions=True)
    fc_extended = FastChem(extended_logk=True)
    fc_extended_w_ions = FastChem(extended_logk=True, with_ions=True)

    assert len(fc_default.gases) > 0
    assert len(fc_w_ions.gases) > 0
    assert len(fc_extended.gases) > 0
    assert fc_extended_w_ions.gases == fc_extended.gases


def test_fastchem_ratio_elements():
    """Test Fastchem with ratio elements"""
    fc = FastChem(ratio_elements=["C", "N"], ratios_to_o=[0.1, 0.2])

    assert fc.ratio_elem_abund["C"] == 0.1
    assert fc.ratio_elem_abund["N"] == 0.2


def test_fastchem_retrieval_parameter_metallicity():
    """Test Fastchem with retrieval parameters."""
    fc = FastChem(metallicity=0.1)

    assert fc.metallicity == 0.1
    fc.metallicity = 0.2
    assert fc.metallicity == 0.2

    fitparams = fc.fitting_parameters()

    assert fitparams["metallicity"][2]() == 0.2

    fitparams["metallicity"][3](0.3)

    assert fc.metallicity == 0.3


def test_fastchem_retrival_parameter_o_ratios():
    """Test Fastchem with retrieval parameters."""
    import random

    fc = FastChem(ratio_elements=["C", "N"], ratios_to_o=[0.1, 0.2])

    # Test general retrieval
    fitparams = fc.fitting_parameters()
    assert "C_O_ratio" in fitparams
    assert "N_O_ratio" in fitparams
    assert "S_O_ratio" in fitparams

    assert fitparams["C_O_ratio"][2]() == 0.1
    assert fitparams["N_O_ratio"][2]() == 0.2

    fitparams["C_O_ratio"][3](0.3)

    assert fc.ratio_elem_abund["C"] == 0.3
    assert fitparams["C_O_ratio"][2]() == 0.3
    assert fc.ratio_elem_abund["N"] == 0.2

    elem_read = {}

    # Test all elements have correct values and do not override each other
    for s in fc.ratio_elem_abund:
        fitparams[f"{s}_O_ratio"][3](random.uniform(0.01, 0.99))  # noqa:S311
        elem_read[s] = fitparams[f"{s}_O_ratio"][2]()

    for s in fc.ratio_elem_abund:
        assert fc.ratio_elem_abund[s] == elem_read[s]


def test_fastchem_compute():
    """Test fastchem computation."""
    import numpy as np

    temperature = np.linspace(2000, 1000, 10)
    pressure = np.linspace(1e6, 1e-4, 10)

    fc = FastChem(selected_elements=["H", "He", "O", "C", "N"])

    fc.initialize_chemistry(10, temperature, pressure)

    assert len(fc.gases) > 0
    assert fc.mixProfile.shape == (len(fc.gases), 10)
    assert fc.mu_profile.shape == (10,)

    amu = fc.muProfile * 6.022e26

    # Near H2 amu
    np.testing.assert_almost_equal(amu, 2.3, decimal=1)

    np.testing.assert_array_equal(fc.muProfile, fc.mu_profile)

    np.testing.assert_almost_equal(fc.mixProfile.sum(axis=0), 1.0)


def test_taurex_detect():
    """Test whether taurex can detect fastchem."""
    from taurex.parameter.classfactory import ClassFactory

    cf = ClassFactory()

    fc = cf.find_klass_from_keyword("fastchem")

    assert fc is not None

    fc = cf.find_klass_from_keyword("fastchem3")

    assert fc is not None
