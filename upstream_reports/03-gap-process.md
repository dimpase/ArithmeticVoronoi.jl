`PERFCOMP_general_query` constructs a shell command from unquoted paths and does not check the backend's exit status before reading its output. I reproduced both behaviours with the original query function and a small backend stub.

Tested at `cfcd79dbddf67a494a1c71b776174b190da5209a` (currently `master`), with GAP 4.16.1 on Linux. No C++ build is needed for the reproducer below.

### Observed behaviour

| Backend fixture | Expected | Observed |
| --- | --- | --- |
| Ordinary path, writes `return [123];`, exits 0 | Return `[123]` | Returns `[123]` (control) |
| Same fixture in a directory containing spaces | Return `[123]` | Shell attempts the path truncated at the space; output cannot be read |
| Ordinary path, writes the same readable output, exits 7 | Report child failure with status and diagnostics | Returns `[123]` as success |

A further fixture that exits 7 without output reaches an output-read error instead of reporting the backend status/stderr. In this GAP version that missing-file path produces a generic `CallFuncList` method error after `ReadAsFunction` fails; the original child diagnostic remains in the temporary error file.

### Source and suggested change

[access_points.g, lines 2452–2460](https://github.com/MathieuDutSik/polyhedral_common/blob/cfcd79dbddf67a494a1c71b776174b190da5209a/CI_tests/access_points.g#L2452-L2460) builds `binary + " " + FileN + " 2> " + FileE`, calls `Exec`, then calls `ReadAsFunction(FileO)()` without validating the child status.

Please use argument-based process execution (e.g. GAP's `Process` API or equivalent), check the exit status before reading output, and preserve/include stderr and the command details when the child fails. A readable file alone is insufficient evidence that the computation succeeded. The same treatment should cover input/error paths containing spaces. This request concerns process handling; changing the GAP serialization format is a separate matter.

### Standalone reproducer

Save the following as `reproduce_gap.py` and run:

```sh
python3 reproduce_gap.py /path/to/polyhedral_common
```

The script loads the query function unchanged from the checkout. It mocks only executable discovery, temporary-directory selection and namelist setup, so the failure is isolated from C++ arithmetic and the discovery helper. The stub reads the actual query output filename, writes valid GAP output and exits with the configured status. Diagnostic files are retained.

```python
"""Run from a polyhedral_common checkout: python3 reproduce_gap.py [checkout]."""
from pathlib import Path
import json
import shutil
import subprocess
import sys
import tempfile

repo = Path(sys.argv[1] if len(sys.argv) > 1 else '.').resolve()
gap = shutil.which('gap')
assert gap, 'GAP is required'
root = Path(tempfile.mkdtemp(prefix='pq_'))
source = (repo / 'CI_tests/access_points.g').read_text()
start = source.index('PERFCOMP_general_query:=function(arg)')
end = source.index('PERFCOMP_group_generators:=function(arg)', start)
query = root / 'original_query.g'
query.write_text(source[start:end])  # unmodified upstream query function

for label, directory, code in [('control', 'bin', 0),
                               ('spaces', 'bin with spaces', 0),
                               ('exit7', 'bin', 7)]:
    work = root / label
    (work / directory).mkdir(parents=True)
    scratch = Path(tempfile.mkdtemp(prefix='pqio_'))
    backend = work / directory / 'backend.py'
    backend.write_text('#!' + sys.executable + '\n' +
        'import pathlib,re,sys\n' +
        'text=pathlib.Path(sys.argv[1]).read_text()\n' +
        'out=re.search(r\'FileListNumberOrbit = "([^\"]+)"\',text).group(1)\n' +
        'pathlib.Path(out).write_text("return [123];\\n")\n' +
        'print("deliberate child diagnostic",file=sys.stderr)\n' +
        f'sys.exit({code})\n')
    backend.chmod(0o700)
    script = work / 'probe.g'
    script.write_text(f'AV_BACKEND := {json.dumps(str(backend))};;\n' +
        f'AV_TMP := {json.dumps(str(scratch))};;\n' +
        'GetBinaryFilename := function(name) return AV_BACKEND; end;\n' +
        '__PERFCOMP_Write_t_space := function(arg) end;\n' +
        'extract_runtime_from_log := function(path) return "unused"; end;\n' +
        'MakeReadWriteGlobal("DirectoryTemporary");\n' +
        'DirectoryTemporary := function() return Directory(AV_TMP); end;\n' +
        f'Read({json.dumps(str(query))});\n' +
        'Print("GAP_VERSION=", GAPInfo.Version, "\\n");\n' +
        'x := PERFCOMP_general_query(rec(), "FileListNumberOrbit");;\n' +
        'Print("QUERY_RETURNED=", x, "\\n");\nQUIT_GAP(0);\n')
    result = subprocess.run([gap, '-q', '-b', '-r', '-A', '-T', '-o', '64m',
                             '-K', '2g', str(script)], capture_output=True,
                            text=True, timeout=30)
    returned = 'QUERY_RETURNED=[ 123 ]' in result.stdout
    print(label, 'intended child exit:', code, 'query returned value:', returned)
    print(result.stdout.strip())
    print(result.stderr.strip())
    child_error = scratch / 'PerfComp.err'
    if child_error.is_file():
        print('Child stderr:', child_error.read_text().strip())
    assert returned == (label != 'spaces')  # observed behaviour at the reported commit

print('Diagnostic files retained at', root)
```

