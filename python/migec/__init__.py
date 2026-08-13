"""migec — UMI barcode extraction, correction and consensus assembly.

The native core is ``migec._core``. It is imported eagerly and on purpose: there is no pure-Python
fallback, because a fallback would make a failed C++ build look like a successful install and only
surface forty minutes into a run.
"""

__version__ = "2.0.0a1"  # keep in sync with pyproject.toml and include/migec/version.hpp

from migec import _core
from migec._core import (
    MIG_FORMAT_VERSION,
    MigecError,
    MigFile,
    MigRecord,
    bucket_of,
    count_fastq,
    pack_barcode,
    reverse_complement,
    unpack_barcode,
)

if _core.__version__ != __version__:
    raise ImportError(
        f"migec version mismatch: python package {__version__} but extension "
        f"{_core.__version__}. Rebuild with `bash setup.sh`."
    )

__all__ = [
    "MIG_FORMAT_VERSION",
    "MigFile",
    "MigRecord",
    "MigecError",
    "__version__",
    "bucket_of",
    "count_fastq",
    "pack_barcode",
    "reverse_complement",
    "unpack_barcode",
]
