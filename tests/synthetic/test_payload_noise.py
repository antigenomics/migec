"""A first-read payload draft is a noisy observation, not an error-free template."""
import csv
import random

import pytest

pytest.importorskip("migec._core", reason="the C++ extension is not built")

from migec.refine import run


@pytest.mark.parametrize("error_quality, corrected", [("+", True), ("I", False)])
def test_payload_compatibility_distinguishes_noise_from_template_evidence(
    tmp_path, error_quality, corrected
):
    rng = random.Random(1911)
    parent, child = 'AAAACCCCGGGG', 'AAAACCCCGGGT'
    template = ''.join(rng.choice('ACGT') for _ in range(32))
    damaged = ''.join(('C' if base == 'A' else 'A') if i < 5 else base
                      for i, base in enumerate(template))
    source = tmp_path / 'tagged.fq'
    with source.open('w') as out:
        for i in range(12):
            out.write(f'@parent{i} RX:Z:{parent}\tQX:Z:IIIIIIIIIIII\n'
                      f'{template}\n+\n' + 'I' * 32 + '\n')
        out.write(f'@child RX:Z:{child}\tQX:Z:IIIIIIIIIIII\n'
                  f'{damaged}\n+\n' + error_quality * 5 + 'I' * 27 + '\n')
        used = {parent, child}
        while len(used) < 3002:
            barcode = ''.join(rng.choice('ACGT') for _ in range(12))
            if barcode in used:
                continue
            used.add(barcode)
            payload = ''.join(rng.choice('ACGT') for _ in range(32))
            out.write(f'@background{len(used)} RX:Z:{barcode}\tQX:Z:IIIIIIIIIIII\n'
                      f'{payload}\n+\n' + 'I' * 32 + '\n')
    run(source, tmp_path / 'out', sample_id='neutral', threads=2)
    with (tmp_path / 'out/neutral.barcodes.tsv').open() as handle:
        table = {row['umi']: row for row in csv.DictReader(handle, delimiter='\t')}
    assert table[child]['parent'] == (parent if corrected else '.')
