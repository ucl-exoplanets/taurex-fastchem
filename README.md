# TauREx-FastChem plugin #

A Python wrapper built using the [TauREx](https://github.com/ucl-exoplanets/taurex3) is available.
The wrapper also installs all available datafiles included with FastChem

*Note: This plugin does not support condensation yet but will do in future versions.*

## Installation


You can install one of the prebuilt binary wheels for Windows, macOS and manylinux through pip:
```bash
pip install taurex-fastchem
```

### Installing from source


To install from source can be done so like this:
```bash
git clone https://github.com/ucl-exoplanets/taurex-fastchem
cd taurex-fastchem
pip install .
```

## Running in TauREx

Once installed you can select the chemical model through the **chemistry_type** keyword under
Chemistry.
```
[Chemistry]
chemistry_type = fastchem
metallicity = 1.0
selected_elements = H, He, C, N, O, Ti, V, S, K
ratio_elements = C, N, Ti
ratios_to_O = 0.5,0.001, 1e-4
with_ions = True

[Fitting]
Ti_O_ratio:fit = True
Ti_O_ratio:prior = "LogUniform(bounds=(-6,2))"
S_O_ratio:fit = True
S_O_ratio:prior = "LogUniform(bounds=(-6,2))"
metallicity:fit = True
metallicity:prior = "LogUniform(bounds=(-6,2))"
```

### Input arguments:

These arguments apply to both the TauREx input file and python interface.

|Argument| Description| Type| Default | Required |
---------|------------|-----|---------|----------|
H_He_ratio| He/H ratio | float | 0.083 | |
selected_elements| List of elements to include in model | list of string | All elements in FastChem | |
ratio_elements| List of elements to set the ratio | list of string | | |
ratios_to_O| ratio of each 'ratio_element' relative to oxygen | array | | |
elements_abundance_file| Path to file that defines initial abundances (in dex) | string | Builtin (solar) | |
metallicity| Metallicity relative to initial abundance | float | 1.0 | |
elements_datafile| Path to file containing elements and their masses | string | Built-in (chemical_abundances.dat) | |
species_datafile| Path to file containing species and thermochemical data | string | Built-in (logK.dat) | |
chem_accuracy| | | | |
with_ions| Include ions | bool | False | |
pressure_accuracy| | | | |
newton_error| | | | |
max_chem_iter| | | | |
max_press_iter| | | | |
max_nedler_iter| | | | |
longdouble| Unused, for compatibility | bool | False | |

### Retrieval Parameters:

|Fitting Parameter| Description| 
---------|------------|
metallicity|Metallicity relative to solar|

The wrapper will generate oxygen retrieval parameters for all metallic elements within the
chemical model. If Ti is present (either by default or specifing in **selected_elements**)
then a **Ti_O_ratio** retrieval parameter will be available.
Using the default **selected_parameters** will give access to:

|Fitting Parameter| Description| 
---------|------------|
Al_O_ratio | Al/O ratio | 
Ar_O_ratio | Ar/O ratio | 
C_O_ratio | C/O ratio | 
Ca_O_ratio | Ca/O ratio | 
Cl_O_ratio | Cl/O ratio | 
Co_O_ratio | Co/O ratio |             
Cr_O_ratio | Cr/O ratio | 
Cu_O_ratio | Cu/O ratio | 
F_O_ratio | F/O ratio | 
Fe_O_ratio | Fe/O ratio | 
Ge_O_ratio | Ge/O ratio | 
K_O_ratio | K/O ratio | 
Mg_O_ratio | Mg/O ratio | 
Mn_O_ratio | Mn/O ratio | 
N_O_ratio | N/O ratio | 
Na_O_ratio | Na/O ratio |
Ne_O_ratio | Ne/O ratio | 
Ni_O_ratio | Ni/O ratio |  
P_O_ratio | P/O ratio | 
S_O_ratio | S/O ratio | 
Si_O_ratio | Si/O ratio |
Ti_O_ratio | Ti/O ratio | 
V_O_ratio | V/O ratio |
Zn_O_ratio | Zn/O ratio |


## Running in Python

You can import the chemistry scheme in Python pretty easily

```python
>>> from pyfastchem import FastChem
>>> fc = FastChem(selected_elements=['H','He','C','O','N','K','e-'], 
                  with_ions=True, metallicity=1.0)
```
You can either pass it into a TauREx forward model like so:
```python
>>> tm = TransmissionModel(chemistry=fc)
```
Or use it independently to compute volume mixing ratios by passing in
temperature and pressure ( Pascal ) arrays:
```python
>>> nlayers = 10
>>> temperature = np.linspace(300,100,nlayers)
>>> pressure = np.logspace(5,-3, nlayers) # Pa
>>> fc.initialize_chemistry(nlayers,temperature,pressure)
>>> fc.gases
['H', 'He', 'O', 'C', 'K', 'N', 'e-', ..., 'O+', 'O-', 'O2+', 'O2-']
>>> fc.mixProfile
array([[3.87435866e-036, 9.95149979e-039, 7.62616463e-042,
        1.23490910e-045, 2.58839801e-050, 3.41640407e-056,
        9.40930967e-064, 9.08433703e-074, 1.41255491e-087,
        1.38065040e-167],
        ...,
       [1.42400626e-001, 1.42400626e-001, 1.42400626e-001,
        1.42400626e-001, 1.42400626e-001, 1.42400791e-001,
        1.42398731e-001, 1.42398284e-001, 1.42367067e-001,
        9.96186945e-001]])
```

