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
