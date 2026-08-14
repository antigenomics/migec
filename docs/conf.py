import os
import pathlib
import sys
import tomllib

sys.path.insert(0, os.path.abspath("../python"))

project = "migec"
copyright = "2026, Mikhail Shugay"
author = "Mikhail Shugay"
# Read, never copied. This was the FOURTH place the version was written by hand and the only one
# nothing checked, which is how every published page said 2.0.0.dev0 for three releases. Importing
# the package instead would not work: `migec._core` is mocked below, so the package's own
# version-agreement assertion fires against a Mock. The file is right there; read it.
release = tomllib.loads(
    (pathlib.Path(__file__).parent.parent / "pyproject.toml").read_text()
)["project"]["version"]

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
    # Two levels open in the sidebar, so every command and every method page is visible without a
    # click. The header carries the seven SECTION names only -- a flat toctree put all twenty page
    # titles up there, and "Fragmented libraries -- why 10x needs a different consensus" is not a
    # navbar link, it is a sentence.
    "show_nav_level": 2,
    "navigation_depth": 2,
    "header_links_before_dropdown": 7,
    "collapse_navigation": False,
    "navbar_align": "left",
}
exclude_patterns = ["_build"]
