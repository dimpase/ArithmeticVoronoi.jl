"""Small, serial runtime reproductions for the three upstream reports.

The boundary-check finding is source-level; this script does not manufacture
an incorrect mathematical boundary or claim an observed nonzero d^2.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import resource
import shutil
import subprocess
import sys
import tempfile
import time

HERE=Path(__file__).resolve().parent
DEFAULT=Path('/tmp/hilbert-ford/arithmetic_voronoi_backend')


def limit():
    resource.setrlimit(resource.RLIMIT_AS,(10*1024**3,10*1024**3))
    resource.setrlimit(resource.RLIMIT_CORE,(0,0))


def execute(argv,folder):
    start=time.monotonic()
    with (folder/'stdout.log').open('w') as out,(folder/'stderr.log').open('w') as err:
        result=subprocess.run(argv,cwd=folder,stdout=out,stderr=err,
                              timeout=45,preexec_fn=limit)
    report=dict(argv=[str(a) for a in argv],exitcode=result.returncode,elapsed_seconds=time.monotonic()-start,
                virtual_memory_limit_bytes=10*1024**3,deadline_seconds=45)
    (folder/'process.json').write_text(json.dumps(report,indent=2)+'\n')
    return report


def official_probe(binary,root):
    public_issue=(HERE/'01-compute-complex.md').read_text()
    input_template=re.search(r'```text\n(.*?)\n```',public_issue,re.S).group(1)+'\n'
    cases={}
    for flag in ['F','T']:
        folder=root/('compute_complex_'+flag);folder.mkdir()
        text=re.sub(r'FileCells = "[^"]*"','FileCells = "cells.g"',input_template)
        text=re.sub(r'FileListNumberOrbit = "[^"]*"','FileListNumberOrbit = "orbit_counts.g"',text)
        text=text.replace('ComputeComplex = F','ComputeComplex = '+flag)
        source=folder/'input.nml';source.write_text(text)
        report=execute([str(binary),str(source)],folder)
        assert report['exitcode']==0
        report['orbit_counts_text']=(folder/'orbit_counts.g').read_text().strip()
        report['executable_sha256']=hashlib.sha256(binary.read_bytes()).hexdigest()
        assert report['orbit_counts_text']=='return [2,4,5,4,1,0];'
        cases[flag]=report
    return cases


def gap_probe(source,gap,root):
    upstream=(source/'CI_tests/access_points.g').read_text()
    begin=upstream.index('PERFCOMP_general_query:=function(arg)')
    end=upstream.index('PERFCOMP_group_generators:=function(arg)',begin)
    function_file=root/'original_general_query.g'
    function_file.write_text(upstream[begin:end])
    cases={}
    for name,with_spaces,exitcode in [('ordinary_success',False,0),('path_with_spaces',True,0),('nonzero_with_output',False,7),('nonzero_without_output',False,7)]:
        folder=root/name;folder.mkdir()
        bin_dir=folder/('bin with spaces' if with_spaces else 'bin');bin_dir.mkdir()
        backend=bin_dir/'probe_backend.py'
        body='''import pathlib,re,sys
text=pathlib.Path(sys.argv[1]).read_text()
output=re.search(r'FileListNumberOrbit = "([^"]+)"',text).group(1)
print('deliberate backend diagnostic',file=sys.stderr)
'''
        if name!='nonzero_without_output':
            body+='pathlib.Path(output).write_text("return [123];\\n")\n'
        body+=f'sys.exit({exitcode})\n'
        backend.write_text('#!'+sys.executable+'\n'+body);backend.chmod(0o700)
        # Keep GAP's generated namelist paths short, as with DirectoryTemporary().
        scratch=Path(tempfile.mkdtemp(prefix='avq_',dir='/tmp'))
        script=folder/'probe.g'
        script.write_text(f'AV_BACKEND := {json.dumps(str(backend))};;\nAV_TMP := {json.dumps(str(scratch))};;\n'+
            '''# Load the original query function; mock only discovery and namelist setup.
GetBinaryFilename := function(name) return AV_BACKEND; end;
__PERFCOMP_Write_t_space := function(arg) end;
extract_runtime_from_log := function(path) return "unused"; end;
MakeReadWriteGlobal("DirectoryTemporary");
DirectoryTemporary := function() return Directory(AV_TMP); end;
'''+
            f'Read({json.dumps(str(function_file))});\n'+
            '''Print("GAP_VERSION=", GAPInfo.Version, "\\n");
answer := PERFCOMP_general_query(rec(), "FileListNumberOrbit");;
Print("QUERY_RETURNED=", answer, "\\n");
QUIT_GAP(0);
''')
        report=execute([gap,'-q','-b','-r','-A','-T','-o','64m','-K','2g',str(script)],folder)
        stdout=(folder/'stdout.log').read_text()
        report['wrapper_returned_value']='QUERY_RETURNED=[ 123 ]' in stdout
        report['intended_backend_exitcode']=exitcode
        report['backend_called']=report['wrapper_returned_value'] or (scratch/'PerfComp.err').is_file() and 'deliberate backend diagnostic' in (scratch/'PerfComp.err').read_text()
        # The success case and the nonzero child with readable output both return normally.
        assert report['wrapper_returned_value']==(name in ['ordinary_success','nonzero_with_output']), (name,stdout,(folder/'stderr.log').read_text())
        if name=='path_with_spaces':
            diagnostic=(scratch/'PerfComp.err').read_text()
            assert any(message in diagnostic for message in ['not found','No such file or directory'])
        if name=='nonzero_without_output':
            assert (scratch/'PerfComp.err').read_text().strip()=='deliberate backend diagnostic'
        shutil.copytree(scratch,folder/'child_files')
        cases[name]=report
    return cases


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--source',type=Path,default=DEFAULT/'polyhedral_common')
    parser.add_argument('--binary',type=Path,default=DEFAULT/'build/bin/PERF_SerialPerfectComputation')
    parser.add_argument('--gap',default=shutil.which('gap'))
    args=parser.parse_args()
    assert args.gap,'GAP is required'
    root=Path(tempfile.mkdtemp(prefix='reproduction_',dir=HERE))
    report=dict(source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=args.source,text=True).strip(),
                official=official_probe(args.binary.resolve(),root),
                gap=gap_probe(args.source.resolve(),args.gap,root),
                boundary_finding='source inspection only; no actual invalid boundary claimed',
                diagnostics=str(root))
    (HERE/'reproduction_summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
