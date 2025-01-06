"""Typing for FastChem objects."""

import typing as t
import numpy.typing as npt
import numpy as np

# struct FastChemInput
# {
#   std::vector<double> temperature;
#   std::vector<double> pressure;

#   bool equilibrium_condensation = false;
#   bool rainout_condensation = false;
# };


class FastChemInput(t.NamedTuple):
    """FastChem input data."""

    temperature: npt.NDArray[np.float64]
    pressure: npt.NDArray[np.float64]
    equilibrium_condensation: bool = False
    rainout_condensation: bool = False


class FastChemOutput(t.NamedTuple):
    """FastChem output data."""

    number_densities: npt.NDArray[np.float64]
    total_element_density: npt.NDArray[np.float64]
    mean_molecular_weight: npt.NDArray[np.float64]
    number_densities_cond: npt.NDArray[np.float64]
    element_cond_degree: npt.NDArray[np.float64]
    element_conserved: npt.NDArray[np.uint32]
    nb_chemistry_iterations: npt.NDArray[np.uint32]
    nb_cond_iterations: npt.NDArray[np.uint32]
    nb_iterations: npt.NDArray[np.uint32]
    fastchem_flag: npt.NDArray[np.uint32]
