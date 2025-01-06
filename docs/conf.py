"""Sphinx configuration."""
project = "TauRex FastChem"
author = "Ahmed F. Al-Refaie"
copyright = "2024, Ahmed F. Al-Refaie"
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx_click",
    "myst_parser",
]
autodoc_typehints = "description"
html_theme = "furo"
