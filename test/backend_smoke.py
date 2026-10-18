"""CTest checks for the native exporter, independent of Julia and saved caches."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile


def namelist(discriminant):
    models = {1: (0, 0), 5: (1, -1), 8: (0, -2),
              12: (0, -3), 13: (1, -3), 40: (0, -10)}
    root_sum, root_product = models[discriminant]
    kind, comm = ('Classic', 'Trivial') if discriminant == 1 else ('RealQuad', 'Use_realimag')
    return f'''&TSPACE
  TypeTspace = "{kind}"
  ClassicDim = 2
  RealImagDim = 2
  RealImagSum = {root_sum}
  RealImagProd = {root_product}
  SuperMatMethod = "NotNeeded"
  ListComm = "{comm}"
  PtGroupMethod = "Trivial"
  FileListSubspaces = "unset"
/
&DATA
  arithmetic = "gmp"
  OnlyWellRounded = T
  ComputeComplex = F
  ComputeBoundary = F
  ComputePreImage = F
  CacheFile = "none"
/
'''


def check(binary, case, folder):
    discriminant = 1 if case in ('classical2', 'rejection') else int(case[1:])
    text = namelist(discriminant)
    if case == 'rejection':
        text = text.replace('ComputeBoundary = F', 'ComputeBoundary = T')
    source, output = folder / 'input.nml', folder / 'output.json'
    source.write_text(text)
    token = 'a' * 64
    result = subprocess.run([binary, str(source), str(output), token],
                            cwd=folder, capture_output=True, text=True, timeout=160)
    if case == 'rejection':
        assert result.returncode == 64, result.stderr
        assert not output.exists(), 'Rejected operation produced an output file'
        errors = [json.loads(line[len('AV_ERROR '):]) for line in result.stderr.splitlines()
                  if line.startswith('AV_ERROR ')]
        assert errors[-1]['code'] == 'unsupported_operation', errors
        return
    assert result.returncode == 0, result.stderr
    data = json.loads(output.read_text())
    assert data['schema'] == 'av_polyhedral_export_v1'
    assert data['request_token'] == token and data['completed'] is True
    assert data['pending_facets'] == []
    representatives = data['representatives']
    expected = {1: 1, 5: 2, 8: 2, 12: 3, 13: 9, 40: 52}[discriminant]
    assert len(representatives) == expected, (len(representatives), expected)
    for representative in representatives:
        assert int(representative['stabilizer_order']) > 0
        assert representative['stabilizer_generators_row']
        assert representative['neighbours']
        for neighbour in representative['neighbours']:
            assert 1 <= neighbour['target'] <= expected
            assert neighbour['facet_vector_indices_zero_based']
            assert neighbour['transporter_row']
    progress = [json.loads(line) for line in
                Path(str(output) + '.progress.jsonl').read_text().splitlines()]
    assert len(progress) == expected and progress[-1]['completed_objects'] == expected
    print(f'{case}: {expected} completed perfect-form orbits')


if __name__ == '__main__':
    binary = str(Path(sys.argv[1]).resolve())
    # Spaces also exercise the argument-based CLI used by the Julia adapter.
    with tempfile.TemporaryDirectory(prefix='av backend smoke ') as temporary:
        check(binary, sys.argv[2], Path(temporary))
