"""Seal the local bundle after validate_fields_delivery.jl and report_fields.py.

Requires the preserved baseline archive and the successful full validation log.
No publication or upload is performed.
"""
from pathlib import Path
import hashlib
import json
import shutil
import tarfile

p=Path(__file__).resolve().parents[1]
r=p/'results';bundle=r/'dataset/v1'
verification=json.loads((r/'fields_dataset_complete_verification.json').read_text())
words=json.loads((r/'fields_words_tests.json').read_text())
assert verification['status']=='passed' and len(verification['examples'])==5
assert words['status']=='passed' and len(words['examples'])==5
log=(r/'fields_delivery_validation.log').read_text()
assert 'ERROR:' not in log and 'Test Failed' not in log
counts=[]
for line in log.splitlines():
    if '|' in line:
        fields=line.split('|',1)[1].split()
        if len(fields)>=2 and fields[0].isdigit() and fields[0]==fields[1]: counts.append(int(fields[0]))
assert counts==[97,6,128,824,2751],counts
total=sum(counts)
summary={'status':'passed','passing_assertions':total,'baseline_assertions':sum(counts[:4]),
         'extended_assertions':counts[-1],'word_examples':sum(x['certified_word_examples'] for x in words['examples'].values()),
         'complete_dataset_replay':verification,'sanitizers':json.loads((r/'sanitizer_fields/summary.json').read_text())}
(r/'fields_delivery_summary.json').write_text(json.dumps(summary,indent=2)+'\n')

def digest(path):
    with path.open('rb') as f: return hashlib.file_digest(f,'sha256').hexdigest()

data=json.loads((bundle/'dataset.json').read_text())
def add(source,name):
    dst=bundle/name;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source,dst)
    ref={'path':name,'sha256':digest(dst),'bytes':dst.stat().st_size};data['artifacts'][name]=ref;return ref

for field in data['fields']:
    D=field['discriminant'];source=r/f'fields/D{D}/PRESENTATIONS.txt'
    field['artifacts']['PRESENTATIONS.txt']=add(source,f'fields/D{D}/PRESENTATIONS.txt')
    cold=r/f'fields/D{D}/initial_cold_stage_checkpoint.json'
    if cold.exists(): field['artifacts']['initial_cold_stage_checkpoint.json']=add(cold,f'fields/D{D}/initial_cold_stage_checkpoint.json')
    if D!=40: add(r/f'fields/D{D}/KW/comparison.json',f'fields/D{D}/KW/comparison.json')
for source in (p/'src').glob('*.jl'): add(source,'software/src/'+source.name)
for name in ['PRESENTATION_JUSTIFICATION.md','CONSTRUCTIVE_WORDS_STATUS.md','KW_COMPARISON.md','DATASET_SCHEMA.md','DATASET_PROPOSAL.md']:
    add(p/name,'documentation/'+name)
for name in ['FIELDS_WORDS_DATA_REPORT.txt','fields_words_tests.json','fields_delivery_summary.json',
             'fields_dataset_complete_verification.json','fields_delivery_validation.log']:
    add(r/name,('documentation/' if name.endswith('.txt') else 'verification/')+name)
backendrefs=[]
for name in ['backend/export_perfect.cpp','backend/pins.json','backend/NOTICE.txt','backend/UPSTREAM_COPYING',
             'backend/PERMUTALIB_COPYING.LGPL','CMakeLists.txt','CMakePresets.json',
             'scripts/run_sanitizer.py','scripts/validate_sanitizer.jl']:
    backendrefs.append(add(p/name,'provenance/'+name))
for source in sorted((p/'cmake').glob('*.cmake*')):
    backendrefs.append(add(source,'provenance/cmake/'+source.name))
add(r/'build/AV_ExportPerfectOrbits.json','provenance/release_build.json')
add(r/'build/AV_ExportPerfectOrbits_sanitizer.json','provenance/sanitizer_build.json')
for source in (p/'test').glob('*.jl'): add(source,'software/test/'+source.name)
data['provenance']['backend_local_source_and_recipe']=backendrefs
data['provenance']['backend_distribution']='external pinned dependency; no upstream source trees or backend binaries bundled'
data['provenance']['validation_assertions']=total
data['provenance']['complete_verification_receipt']='verification/fields_dataset_complete_verification.json'
(bundle/'dataset.json').write_text(json.dumps(data,indent=2)+'\n')
assert all(digest(bundle/ref['path'])==ref['sha256'] for ref in data['artifacts'].values())

baseline=json.loads((r/'baselines/two_cells_2026_10_17/manifest.json').read_text())
critical=[name for name in baseline['file_sha256'] if name.startswith('results/two_cells/')]
critical += ['src/Validation.jl','Project.toml','Manifest.toml','backend/config.json',
             'backend/bin/AV_ExportPerfectOrbits','BACKEND_AUDIT.md','results/D5.json','results/D8.json']
assert all(digest(p/name)==baseline['file_sha256'][name] for name in critical)
(r/'fields_baseline_preservation.json').write_text(json.dumps({'status':'passed','checked_unchanged':critical,
    'snapshot_archive':'results/baselines/two_cells_2026_10_17/baseline.tar.gz'},indent=2)+'\n')

archive=r/'dataset/arithmetic_groups_dataset_v1.tar.gz'
with tarfile.open(archive,'w:gz') as tar: tar.add(bundle,arcname='arithmetic_groups_dataset_v1')
receipt={'archive':archive.name,'sha256':digest(archive),'bytes':archive.stat().st_size,
         'artifacts':len(data['artifacts']),'dataset_manifest_sha256':digest(bundle/'dataset.json'),
         'publication_status':'local candidate, no DOI, no submission or contact'}
(archive.parent/'archive_manifest.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(total,'passing assertions;',len(data['artifacts']),'dataset artifacts;',archive.stat().st_size,'archive bytes')
