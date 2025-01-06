"""TauRex FastChem."""

import pyfastchem
from .types import FastChemInput, FastChemOutput
import typing as t
import os
from taurex.chemistry import AutoChemistry
from taurex.core import fitparam
import numpy as np
import importlib.resources as ires
from .util import _create_pyfastchem_, fix_fastchem_species

DEFAULT_LOGK_WO_IONS = ires.files("taurex_fastchem") / "data" / "logK_wo_ions.dat"
DEFAULT_LOGK = ires.files("taurex_fastchem") / "data" / "logK.dat"
DEFAULT_LOGK_EXTENDED = ires.files("taurex_fastchem") / "data" / "logK_extended.dat"

CONST_K = 1.3806504e-16


class FastChem(AutoChemistry):
    """FastChem plugin wrapper for TauREx."""

    def __init__(
        self,
        h_he_ratio: t.Optional[float] = None,
        selected_elements: t.Optional[t.List[str]] = None,
        elements_abundance_file: t.Optional[os.PathLike] = None,
        ratio_elements: t.Optional[t.List[str]] = None,
        ratios_to_o: t.Optional[t.List[float]] = None,
        metallicity: t.Optional[float] = 1.0,
        elements_datafile: t.Optional[os.PathLike] = None,
        species_datafile: t.Optional[os.PathLike] = None,
        chem_accuracy: t.Optional[float] = 1e-4,
        with_ions: t.Optional[bool] = False,
        extended_logk: t.Optional[bool] = False,
        pressure_accuracy: t.Optional[float] = 1e-4,
        newton_error: t.Optional[float] = 1e-4,
        max_chem_iter: t.Optional[int] = 300,
        max_press_iter: t.Optional[int] = 100,
        max_nedler_iter: t.Optional[int] = 100,
        longdouble: t.Optional[bool] = False,
    ):
        """Initialise the FastChem plugin."""
        super().__init__(f"FastChem_{3.0}")

        self._gases = None
        self._mixprofile = None

        self._fc_input_data: FastChemInput = pyfastchem.FastChemInput()
        self._fc_output_data: FastChemOutput = pyfastchem.FastChemOutput()

        self._elements_datafile = elements_datafile
        self._species_datafile = species_datafile
        self._chem_accr = chem_accuracy
        self._press_accuracy = pressure_accuracy
        self._newton_error = newton_error
        self._max_chem_iter = max_chem_iter
        self._max_press_iter = max_press_iter
        self._max_nedler_iter = max_nedler_iter
        self._with_ions = with_ions
        self._metallicity = metallicity
        self._h_he_ratio = h_he_ratio

        self._log_k_file = self._species_datafile

        if self._log_k_file is None:
            if extended_logk:
                self._log_k_file = DEFAULT_LOGK_EXTENDED
            elif with_ions:
                self._log_k_file = DEFAULT_LOGK
            else:
                self._log_k_file = DEFAULT_LOGK_WO_IONS

        self.selected_elements = selected_elements
        self.element_abundance_file = elements_abundance_file

        self._ratio_elements = ratio_elements
        self._ratios_to_o = ratios_to_o

        self._fc = _create_pyfastchem_(
            self._log_k_file, self.selected_elements, self.element_abundance_file
        )

        element_count = self._fc.getElementNumber()

        elements = [self._fc.getElementSymbol(i) for i in range(element_count)]
        abundances = [self._fc.getElementAbundance(i) for i in range(element_count)]

        self.elem_abundance = dict(zip(elements, abundances))

        self._h_he_ratio = h_he_ratio or self.elem_abundance["He"]

        self.elem_abundance["He"] = self._h_he_ratio

        self.info("Current elements: %s", self.elem_abundance)

        # ratios_elements = ratio_elements or []
        # ratio_abundances = ratios_to_o or [
        #     self.elem_abundance[x] / self.elem_abundance["O"] for x in ratios_elements
        # ]

        all_ratio_elements = [
            x
            for x in self.elem_abundance.keys()
            if x
            not in (
                "H",
                "He",
                "O",
                "el",
            )
        ]
        ratio_abundances = [
            self.elem_abundance[x] / self.elem_abundance["O"]
            for x in all_ratio_elements
        ]

        self.ratio_elem_abund = dict(zip(all_ratio_elements, ratio_abundances))

        if ratio_elements is not None:
            for elem, abund in zip(ratio_elements, ratios_to_o):
                self.ratio_elem_abund[elem] = abund

        self.info("Current ratios: %s", self.ratio_elem_abund)

        num_species = self._fc.getGasSpeciesNumber()
        self.gas_species = [self._fc.getGasSpeciesSymbol(i) for i in range(num_species)]

        self.gas_species = fix_fastchem_species(self.availableActive, self.gas_species)
        self._mixprofile = None
        self.add_ratio_parameters()

    def add_ratio_parameters(self):
        """Add ratio to retrievsal fitting parameters."""
        for elem in self.ratio_elem_abund.keys():
            param_name = f"{elem}_O_ratio"
            param_latex = f"[{elem}/O]"

            def read_mol(self, elem=elem):
                return self.ratio_elem_abund[elem]

            def write_mol(self, value, elem=elem):
                self.ratio_elem_abund[elem] = value

            read_mol.__doc__ = f"{elem}/O ratio"

            bounds = [1e-12, 0.1]

            self.add_fittable_param(
                param_name, param_latex, read_mol, write_mol, "log", False, bounds
            )

    @property
    def gases(self):
        """Return the gas species symbols."""
        return self.gas_species

    @property
    def mixProfile(self):  # noqa: N802
        """Return the computed mixing profile."""
        return self._mixprofile

    @fitparam(
        param_name="metallicity",
        param_latex="Z",
        default_bounds=[0.2, 2.0],
        default_fit=False,
    )
    def metallicity(self):
        """Metallicity parameter."""
        return self._metallicity

    @metallicity.setter
    def metallicity(self, value):
        """Set the metallicity parameter."""
        self._metallicity = value

    def update_abundances(self):
        """Update the abundances for fastchem."""
        # Compute new O
        new_o = self.elem_abundance["O"] * self._metallicity
        new_he = self._h_he_ratio

        new_elem_abund = [1.0, new_he, new_o]  # Hydrogen
        for elem, abund in self.elem_abundance.items():
            if elem in ("H", "He", "O"):
                continue
            new_elem_abund.append(self.ratio_elem_abund[elem] * new_o)

        self._fc.setElementAbundances(new_elem_abund)

    def initialize_chemistry(
        self, nlayers, temperature_profile, pressure_profile, altitude_profile=None
    ):
        """Initialize the chemistry."""
        self.update_abundances()

        density = pressure_profile * 10 / (CONST_K * temperature_profile)

        self._fc_input_data.temperature = temperature_profile
        # Pa to BAr
        self._fc_input_data.pressure = pressure_profile / 1e5

        self._fc_input_data.equilibrium_condensation = False
        self._fc_input_data.rainout_condensation = False

        self._fc.calcDensities(self._fc_input_data, self._fc_output_data)

        number_densities = self._fc_output_data.number_densities

        self._mixprofile = np.array(number_densities).T / density

        sum_mix = self._mixprofile.sum(axis=0)

        self._mixprofile /= sum_mix

        mean_molecular_weight = np.array(self._fc_output_data.mean_molecular_weight)
        # Conver it to kg from g/mol

        self.mu_profile = mean_molecular_weight * 1e-3 / 6.02214076e23

    @classmethod
    def input_keywords(cls):
        return ["fastchem", "fastchem3"]

    BIBTEX_ENTRIES = [
        """
        @article{fastchem,
            author = {Stock, Joachim W and Kitzmann, Daniel and Patzer, A Beate C and Sedlmayr, Erwin},
            title = "{FastChem: A computer program for efficient complex chemical equilibrium calculations in the neutral/ionized gas phase with applications to stellar and planetary atmospheres }",
            journal = {Monthly Notices of the Royal Astronomical Society},
            volume = {479},
            number = {1},
            pages = {865-874},
            year = {2018},
            month = {06},
            abstract = "{For the calculation of complex neutral/ionized gas-phase chemical equilibria, we present a semi-analytical, versatile, and efficient computer program, called FastChem. The applied method is based on the solution of a system of coupled non-linear (and linear) algebraic equations, namely the law of mass action and the element conservation equations including charge balance, in many variables. Specifically, the system of equations is decomposed into a set of coupled nonlinear equations in one variable each, which are solved analytically whenever feasible to reduce computation time. Notably, the electron density is determined by using the method of Nelder and Mead at low temperatures. The program is written in object-oriented C++ which makes it easy to couple the code with other programs, although a stand-alone version is provided. FastChem can be used in parallel or sequentially and is available under the GNU General Public License version 3 at https://github.com/exoclime/FastChem together with several sample applications. The code has been successfully validated against previous studies and its convergence behaviour has been tested even for extreme physical parameter ranges down to \\$100\\,\\mathrm\\{K\\}\\$ and up to \\$1000\\,\\mathrm\\{bar\\}\\$. FastChem converges stable and robust in even most demanding chemical situations, which posed sometimes extreme challenges for previous algorithms.}",
            issn = {0035-8711},
            doi = {10.1093/mnras/sty1531},
            url = {https://doi.org/10.1093/mnras/sty1531},
            eprint = {https://academic.oup.com/mnras/article-pdf/479/1/865/25126582/sty1531.pdf},
        }
        """
    ]
