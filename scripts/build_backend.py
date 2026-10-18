"""Bounded, serial build of one pinned upstream target; retain commands and logs."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import resource
import re
import shlex
import shutil
import subprocess
import time

PROJECT=Path(__file__).resolve().parents[1]
BASE=Path('/tmp/hilbert-ford/arithmetic_voronoi_backend')


def capture(argv,cwd=None):
    return subprocess.check_output(argv,cwd=cwd,text=True,stderr=subprocess.STDOUT).strip()


def limit():
    resource.setrlimit(resource.RLIMIT_AS,(10*1024**3,10*1024**3))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--target',default='PERF_SerialPerfectComputation')
    parser.add_argument('--build',type=Path,default=BASE/'build')
    parser.add_argument('--source',type=Path,default=BASE/'polyhedral_common')
    parser.add_argument('--sanitizer',action='store_true')
    parser.add_argument('--configure-default',action='store_true')
    parser.add_argument('--configure-fields',action='store_true',help='install the extended-field executable/config without replacing the accepted baseline')
    args=parser.parse_args()
    assert not (args.configure_default and args.configure_fields),'select only one installed configuration'
    assert args.target in ['PERF_SerialPerfectComputation','AV_ExportPerfectOrbits']
    output=PROJECT/'results/build';output.mkdir(parents=True,exist_ok=True)
    label=args.target+('_sanitizer' if args.sanitizer else '')
    binary=args.build/'bin'/label
    compiler='c++'
    if args.sanitizer:
        assert args.target=='AV_ExportPerfectOrbits'
        capture(['cmake','-S',str(args.source),'-B',str(args.build),'-DCMAKE_EXPORT_COMPILE_COMMANDS=ON'])
        commands=json.loads((args.build/'compile_commands.json').read_text())
        entry=next(c for c in commands if Path(c['file']).resolve()==(PROJECT/'backend/export_perfect.cpp').resolve())
        words=entry.get('arguments') or shlex.split(entry['command'])
        includes=[];defines=[];index=0
        while index<len(words):
            word=words[index]
            if word.startswith('-I'):includes.append(word)
            elif word=='-isystem':includes.extend(words[index:index+2]);index+=1
            elif word.startswith('-D') and word!='-DNDEBUG':defines.append(word)
            index+=1
        compiler=shutil.which('clang++')
        assert compiler,'clang++ is needed for this sanitizer build'
        command=[compiler,'-std=c++20','-O1','-gline-tables-only','-UNDEBUG','-fno-omit-frame-pointer',
            '-fsanitize=address,undefined','-Wall','-Wextra',*defines,*includes,
            str(PROJECT/'backend/export_perfect.cpp'),str(args.build/'external/nauty/install/lib/libnauty.a'),
            '-lgmpxx','-lgmp','-lboost_serialization','-pthread','-lm','-o',str(binary)]
    else:
        command=['cmake','--build',str(args.build),'--target',args.target,'--parallel','1']
    log=output/(label+'.log')
    if log.exists():
        history=output/'history'/str(time.time_ns());history.mkdir(parents=True)
        shutil.copy2(log,history/log.name)
        previous=output/(label+'.json')
        if previous.exists():shutil.copy2(previous,history/previous.name)
    start=time.monotonic()
    with log.open('w') as stream:
        process=subprocess.Popen(command,stdout=stream,stderr=subprocess.STDOUT,preexec_fn=limit)
        code=process.wait()
    info=dict(command=command,returncode=code,elapsed_seconds=time.monotonic()-start,
        virtual_memory_limit_bytes=10*1024**3,parallel_build_jobs=1,
        source=str(args.source),top_commit=capture(['git','rev-parse','HEAD'],args.source),
        submodules=capture(['git','submodule','status','--recursive'],args.source).splitlines(),
        source_diff=capture(['git','diff'],args.source),
        compiler=capture([compiler,'--version']).splitlines()[0],cmake=capture(['cmake','--version']).splitlines()[0],
        gmp=capture(['pkg-config','--modversion','gmpxx']),eigen=capture(['pkg-config','--modversion','eigen3']),
        cmake_cache=(args.build/'CMakeCache.txt').read_text(),log=str(log),
        binary=str(binary.resolve()),environment={k:os.environ.get(k) for k in ['CXX','CXXFLAGS','CMAKE_BUILD_PARALLEL_LEVEL']})
    info['exporter_source_sha256']=hashlib.sha256((PROJECT/'backend/export_perfect.cpp').read_bytes()).hexdigest()
    info['nauty_commit']=capture(['git','rev-parse','HEAD'],args.build/'external/src/NautyDirectPackage')
    info['nauty_config_log']=(args.build/'external/src/NautyDirectPackage/config.log').read_text()
    info['boost_version_header']=next(line for line in Path('/usr/include/boost/version.hpp').read_text().splitlines() if line.startswith('#define BOOST_LIB_VERSION'))
    info['sanitizers']='address,undefined; C/C++ assertions enabled; -Wall -Wextra' if args.sanitizer else 'SANITY_CHECK in the release build'
    compile_database=args.build/'compile_commands.json'
    if compile_database.exists():
        info['target_compile_commands']=[c for c in json.loads(compile_database.read_text())
            if Path(c['file']).name in [args.target+'.cpp','export_perfect.cpp']]
    if code==0 and binary.is_file():
        info['executable_sha256']=hashlib.sha256(binary.read_bytes()).hexdigest()
    (output/(label+'.json')).write_text(json.dumps(info,indent=2)+'\n')
    if args.configure_default or args.configure_fields:
        assert code==0 and not args.sanitizer,'Only an explicitly selected successful release build can become default'
        destination=PROJECT/'backend/bin'/(args.target+('_fields' if args.configure_fields else ''))
        destination.parent.mkdir(parents=True,exist_ok=True)
        if destination.resolve()!=binary.resolve():shutil.copy2(binary,destination)
        assert hashlib.sha256(destination.read_bytes()).hexdigest()==info['executable_sha256']
        (PROJECT/'backend'/('config_fields.json' if args.configure_fields else 'config.json')).write_text(json.dumps(dict(schema='av_backend_config_v1',
            executable=str(destination.resolve()),executable_sha256=info['executable_sha256'],commit=info['top_commit']),indent=2)+'\n')
    print(json.dumps({k:v for k,v in info.items() if k not in ['cmake_cache','source_diff','submodules','nauty_config_log']},indent=2))
    raise SystemExit(code)


if __name__=='__main__':main()
