#!/usr/bin/env bash
# Development install. Creates .venv with uv, builds the extension in editable mode, and then
# asserts that the EXTENSION imports -- not just the package. Without that assertion a failed C++
# build looks like a successful install, and the failure only surfaces much later.
set -euo pipefail

cd "$(dirname "$0")"

if ! command -v uv >/dev/null 2>&1; then
  echo "uv not found: https://docs.astral.sh/uv/getting-started/installation/" >&2
  exit 1
fi

[ -d .venv ] || uv venv --python 3.12
# shellcheck disable=SC1091
source .venv/bin/activate

# --no-build-isolation keeps the editable rebuild hook working against this venv's pybind11.
uv pip install scikit-build-core "pybind11>=3.0.2,<4" ninja cmake
uv pip install -e ".[test,dev]" --no-build-isolation

python - <<'PY'
import migec
from migec import _core
assert _core.__version__ == migec.__version__, (_core.__version__, migec.__version__)
print(f"migec {migec.__version__}  core {_core.__version__}  mig format v{_core.MIG_FORMAT_VERSION}")
PY

echo
echo "C++ tests:  cmake -S . -B build -DMIGEC_TESTS=ON && cmake --build build -j && ctest --test-dir build"
echo "Python:     python -m pytest tests/unit tests/synthetic -q"
echo "Docs:       sphinx-build -W --keep-going -b html docs docs/_build/html"
# The notebooks are NOT installed here on purpose: each declares its own dependencies in a PEP 723
# header, so `uv run` builds an environment per notebook and they stay runnable for someone who
# never cloned the repo. Running one with this venv's bare python instead fails on a missing
# polars, which looks like a broken notebook rather than a missing extra.
echo "Notebooks:  uv run marimo edit notebooks/platforms.py   (deps come from the notebook header)"
