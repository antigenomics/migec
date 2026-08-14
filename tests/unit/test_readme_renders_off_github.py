"""The README is the PyPI long description, so a repo-relative link is broken on arrival.

PyPI serves the rendered README from its own domain and does not rewrite relative URLs, so
`src="assets/pipeline.svg"` is a broken image on every project page from the moment the wheel
lands -- and it looks perfect in the GitHub preview, which is where anyone would check. Both
images were relative through 2.3.0.
"""

from __future__ import annotations

import re
from pathlib import Path

README = Path(__file__).resolve().parents[2] / "README.md"


def test_every_image_and_link_is_absolute():
    text = README.read_text()
    srcs = re.findall(r'src="([^"]+)"', text)
    assert srcs, "the README lost its images"
    relative = [s for s in srcs if not s.startswith(("http://", "https://"))]
    assert not relative, f"repo-relative image src, broken on PyPI: {relative}"

    # Markdown links too: `](assets/x.svg)` and `](ROADMAP.md)` are the same trap.
    targets = re.findall(r"\]\(([^)]+)\)", text)
    relative_links = [t for t in targets if not t.startswith(("http://", "https://", "#"))]
    assert not relative_links, f"repo-relative link, broken on PyPI: {relative_links}"
