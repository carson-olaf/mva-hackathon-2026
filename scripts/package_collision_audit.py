"""Build a local public-only replay archive from an explicit file allowlist.

No broad repository walk, patient input, upload, or publication is performed.
"""
import argparse
import hashlib
import json
from pathlib import Path
import zipfile

ROOT=Path(__file__).resolve().parents[1]
SCRIPTS=[
    'public_collision_fetch.py','audit_target_collisions.py',
    'audit_collision_engine_sources.py','build_collision_evidence.py',
    'replay_collision_audit.py','independent_collision_check.py',
    'package_collision_audit.py',
]
RESULTS=[
    'summary.json','engine_source_audit.json','live_api_evidence.json',
    'identity_conflicts.csv','current_identity_conflicts.csv',
    'rare_disease_exposure.csv','rare_disease_gene_evidence.csv','named_cases.csv',
    'input_manifest.json','expected_outputs.json','README.md','NOTICE.md',
]


def package(cache, destination):
    manifest=json.loads((ROOT/'results/collision_audit/input_manifest.json').read_text())
    sources=[(ROOT/'scripts'/n,Path('scripts')/n) for n in SCRIPTS]
    sources += [(ROOT/'results/collision_audit'/n,Path('results/collision_audit')/n) for n in RESULTS]
    sources += [(ROOT/'report'/n,Path('report')/n) for n in ['collision_audit.md','collision_maintainer_note.md','target_identity_audit.md']]
    sources += [(ROOT/'tests/test_collision_audit.py',Path('tests/test_collision_audit.py'))]
    for record in manifest['requests']:
        name=record['filename']
        if Path(name).name!=name:
            raise ValueError('Expected plain public cache filename')
        path=cache/name
        if hashlib.sha256(path.read_bytes()).hexdigest()!=record['sha256']:
            raise ValueError(f'Hash mismatch: {name}')
        sources.append((path,Path('work/public-collision-audit')/name))
    destination.parent.mkdir(parents=True,exist_ok=True)
    temporary=destination.with_suffix('.part')
    with zipfile.ZipFile(temporary,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as archive:
        for path,relative in sources:
            archive.write(path,relative.as_posix())
        archive.writestr('work/public-collision-audit/manifest.json',json.dumps(manifest,indent=2)+'\n')
        archive.writestr('README.md',
            '# Public target-identity audit\n\n'
            'Offline replay: `python3 scripts/replay_collision_audit.py`.\n\n'
            'Independent census: `python3 scripts/independent_collision_check.py`.\n\n'
            'See results/collision_audit/README.md, NOTICE.md and report/collision_audit.md.\n'
            'Only public inputs are included. No patient or gated data.\n')
        archive.writestr('pyproject.toml',
            '[project]\nname = "public-target-identity-audit"\nversion = "0.1.0"\n'
            'requires-python = ">=3.12"\ndependencies = ["requests>=2.32,<3"]\n\n'
            '[dependency-groups]\ndev = ["pytest>=9,<10"]\n\n[tool.uv]\npackage = false\n')
    temporary.replace(destination)
    digest=hashlib.sha256(destination.read_bytes()).hexdigest()
    destination.with_suffix('.zip.sha256').write_text(f'{digest}  {destination.name}\n')
    print(json.dumps({'archive':str(destination),'bytes':destination.stat().st_size,
                      'sha256':digest,'allowlisted_source_files':len(sources)}))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache',type=Path,default=ROOT/'work/public-collision-audit')
    parser.add_argument('--output',type=Path,default=ROOT/'work/public-collision-audit-replay.zip')
    args=parser.parse_args()
    package(args.cache,args.output)
