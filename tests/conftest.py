"""Shared skip conditions.

Gating is done with skipif CONSTANTS rather than pytest markers: a marker has to be registered,
documented and remembered, whereas a constant carries its own reason string telling you the exact
command that fixes it.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest


def _core_available() -> bool:
    try:
        import migec._core  # noqa: F401
    except Exception:
        return False
    return True


requires_core = pytest.mark.skipif(
    not _core_available(),
    reason="the native extension is not built -- run `bash setup.sh`",
)

requires_benchmark = pytest.mark.skipif(
    not os.getenv("RUN_BENCHMARK"),
    reason="set RUN_BENCHMARK=1 to run benchmarks",
)

# The real-data tier reads the `ci/` fixtures of `isalgo/umi_data`. They are not vendored here:
# they are LFS objects, and the mirror is where they are written from anyway (SOURCES.md).
UMI_DATA = Path(os.environ.get("UMI_DATA", "~/hf/umi_data")).expanduser()

requires_umi_data = pytest.mark.skipif(
    not (UMI_DATA / "ci").is_dir(),
    reason=(
        f"the CI fixtures are not at {UMI_DATA} -- `git clone "
        f"https://huggingface.co/datasets/isalgo/umi_data {UMI_DATA}` (git-lfs required), or set "
        f"UMI_DATA to a copy"
    ),
)
