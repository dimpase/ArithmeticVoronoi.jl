"""Configure only the pinned, narrow backend. Never advance submodules remotely."""
import argparse
import json
from pathlib import Path
import subprocess

PROJECT=Path(__file__).resolve().parents[1]
BASE=Path('/tmp/hilbert-ford/arithmetic_voronoi_backend')


def run(argv, **kwargs):
    subprocess.run([str(x) for x in argv],check=True,**kwargs)


def output(argv,cwd):
    return subprocess.check_output(argv,cwd=cwd,text=True).strip()


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--source',type=Path,default=BASE/'polyhedral_common')
    parser.add_argument('--build',type=Path,default=BASE/'build')
    args=parser.parse_args()
    source=args.source.resolve();build=args.build.resolve()
    pins=json.loads((PROJECT/'backend/pins.json').read_text())
    if not source.exists():
        source.parent.mkdir(parents=True,exist_ok=True)
        run(['git','clone','--no-checkout',pins['repository'],source])
        run(['git','checkout','--detach',pins['commit']],cwd=source)
    assert output(['git','rev-parse','HEAD'],source)==pins['commit'],'Refusing an unpinned checkout'
    run(['git','submodule','update','--init','--recursive'],cwd=source)
    lines=output(['git','submodule','status','--recursive'],source).splitlines()
    actual={line.strip().split()[1]:line.strip().split()[0] for line in lines}
    assert actual==pins['submodules'],'Submodule commits differ from the accepted build'
    patch=PROJECT/'patches/polyhedral_common.patch'
    if subprocess.run(['git','apply','--check',str(patch)],cwd=source,capture_output=True).returncode==0:
        run(['git','apply',patch],cwd=source)
    else:
        run(['git','apply','--reverse','--check',patch],cwd=source)
    # Verify the complete tracked diff, rather than merely accepting a patch hunk.
    actual_diff=output(['git','diff','--no-ext-diff'],source)
    assert actual_diff.strip()==patch.read_text().strip(),'Unexpected upstream edits; use an isolated checkout'
    run(['cmake','-S',source,'-B',build,'-G','Ninja',
         '-DWITH_MPI=OFF','-DENABLE_FLINT_SUPPORT=OFF','-DENABLE_SCIP_SUPPORT=OFF',
         '-DENABLE_BOOST_TYPES=OFF','-DGMP_POOL=OFF','-DCMAKE_BUILD_TYPE=Release',
         '-DCMAKE_CXX_FLAGS_RELEASE=-O1 -DNDEBUG','-DCMAKE_CXX_FLAGS=-DSANITY_CHECK',
         '-DCMAKE_EXPORT_COMPILE_COMMANDS=ON',f'-DAV_EXPORTER_SOURCE={PROJECT / "backend/export_perfect.cpp"}'])
    print('Configured pinned source at',source,'and serial build at',build)


if __name__=='__main__':main()
