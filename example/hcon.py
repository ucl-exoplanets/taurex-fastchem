from taurex_fastchem import FastChem
import matplotlib.pyplot as plt
import numpy as np

fc = FastChem(selected_elements=["H", "He", "O", "C", "N"])

temperature = np.linspace(3000, 1000, 100)
pressure = np.logspace(6, -2, 100) * 1e5  # Bar to Pascal

fc.initialize_chemistry(100, temperature_profile=temperature, pressure_profile=pressure)

species_to_see = ["H2", "H2O1", "C1H4", "H3N1", "C2H2", "C1O1", "C1O2", "C1H2O1"]

mixprofile = fc.mixProfile
pressure_bar = pressure / 1e5

amu_kg = fc.mu_profile * 6.022e26

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
for i, species in enumerate(species_to_see):
    ax1.plot(mixprofile[fc.gases.index(species)], pressure_bar, label=species)

ax1.set_yscale("log")
ax1.set_xscale("log")
ax1.invert_yaxis()
ax1.set_ylabel("Pressure [bar]")
ax1.set_xlabel("Mixing ratio")

ax1.legend()

ax2.plot(amu_kg, pressure_bar)
ax2.set_yscale("log")
ax2.invert_yaxis()
ax2.set_ylabel("Pressure [bar]")
ax2.set_xlabel("Mean molecular weight [kg]")
plt.show()
