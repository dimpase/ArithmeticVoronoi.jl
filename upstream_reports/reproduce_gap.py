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
