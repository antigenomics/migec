import os
import sys

sys.path.insert(0, os.path.abspath("../python"))

project = "migec"
copyright = "2026, Mikhail Shugay"
author = "Mikhail Shugay"
release = "2.0.0.dev0"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.githubpages",
    "sphinx.ext.mathjax",
]

napoleon_google_docstring = True
napoleon_numpy_docstring = False
autodoc_member_order = "bysource"
autodoc_typehints = "description"
autodoc_default_options = {"members": True, "undoc-members": True, "show-inheritance": True}
# The extension is not built when the docs are: mock it rather than requiring a compiler in the
# docs job.
autodoc_mock_imports = ["typer", "polars", "migec._core"]

html_theme = "pydata_sphinx_theme"
html_theme_options = {
    "github_url": "https://github.com/antigenomics/migec",
    "show_nav_level": 2,
}
exclude_patterns = ["_build"]
