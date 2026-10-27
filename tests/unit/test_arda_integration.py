"""The optional ARDA module installs its actual distribution and uses real flags."""

from pathlib import Path
import re

import pytest

MODULE = Path(__file__).resolve().parents[2] / "integrations/nextflow/modules/downstream/arda"


def test_environment_installs_arda_mapper():
    environment = (MODULE / "environment.yml").read_text()
    assert "- arda-mapper>=2.36.2" in environment
    assert not re.search(r"^\s*- arda\s*$", environment, re.MULTILINE)


def test_module_flags_exist_on_arda_amplicon():
    pytest.importorskip("arda", reason="optional ARDA integration needs arda-mapper")
    import typer.main
    from arda.cli import app

    command = typer.main.get_command(app).commands["amplicon"]
    available = {flag for parameter in command.params for flag in parameter.opts}
    script = (MODULE / "main.nf").read_text().split("script:", 1)[1].split("cat <<-END_VERSIONS", 1)[0]
    assert set(re.findall(r"--[a-z][a-z0-9-]*", script)) <= available


def test_module_does_not_default_to_an_unpublished_image():
    source = (MODULE / "main.nf").read_text()
    assert "getOrDefault('arda_container', null)" in source
