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


@pytest.mark.parametrize('different_mate', [False, True])
@pytest.mark.parametrize('budget, threads', [(0, 1), (4096, 4)])
@pytest.mark.parametrize('prefix_error_quality', ['+', '!', 'I'])
def test_paired_context_protects_distinct_template_after_noisy_prefix(
    tmp_path, different_mate, budget, threads, prefix_error_quality
):
    from migec.checkout import run as checkout_run

    rng = random.Random(281)
    parent, child = 'AAAACCCCGGGG', 'AAAACCCCGGGT'
    adapter = 'CAGTGGTATCAACGCAGAGT'
    template = ''.join(rng.choice('ACGT') for _ in range(32))
    mate = ''.join(rng.choice('ACGT') for _ in range(128))
    damaged = ''.join(('C' if b == 'A' else 'A') if i < 5 else b
                      for i, b in enumerate(template))
    child_mate = ''.join(('C' if b == 'A' else 'A') if different_mate and 40 <= i < 60 else b
                         for i, b in enumerate(mate))
    r1, r2 = tmp_path / 'r1.fq', tmp_path / 'r2.fq'
    with r1.open('w') as a, r2.open('w') as b:
        def write(name, barcode, payload, payload_quality, mate_sequence):
            sequence = barcode + adapter + payload
            quality = 'I' * (len(barcode) + len(adapter)) + payload_quality
            a.write(f'@{name}\n{sequence}\n+\n{quality}\n')
            b.write(f'@{name}\n{mate_sequence}\n+\n' + 'I' * len(mate_sequence) + '\n')

        for i in range(12):
            write(f'parent{i}', parent, template, 'I' * 32, mate)
        write('child', child, damaged, prefix_error_quality * 5 + 'I' * 27, child_mate)
        used = {parent, child}
        while len(used) < 3002:
            barcode = ''.join(rng.choice('ACGT') for _ in range(12))
            if barcode in used:
                continue
            used.add(barcode)
            payload = ''.join(rng.choice('ACGT') for _ in range(32))
            other_mate = ''.join(rng.choice('ACGT') for _ in range(128))
            write(f'background{len(used)}', barcode, payload, 'I' * 32, other_mate)
    sheet = tmp_path / 'barcodes.tsv'
    sheet.write_text('neutral\t^' + 'N' * 12 + adapter.lower() + '\n')
    checkout_run(r1, sheet, tmp_path / 'checkout', reads2=r2, mig=True, threads=2)
    report = run(tmp_path / 'checkout', tmp_path / 'refine', sample_id='neutral',
                 threads=threads, table_budget_bytes=budget)
    assert report['payload_primary_width'] == 32
    assert report['payload_observed_width'] == 160
    assert report['payload_mate1_span'] == 32
    with (tmp_path / 'refine/neutral.barcodes.tsv').open() as handle:
        table = {row['umi']: row for row in csv.DictReader(handle, delimiter='\t')}
    assert table[child]['parent'] == ('.' if different_mate else parent)
