"""Reproduce the public collision audit. Offline by default; no private inputs.

--download explicitly fetches the pinned public requests if they are missing.
Mutable sources must still match the recorded SHA256: drift is an error, not a
successful reproduction. The preserved public cache can be used offline.
"""
import argparse
import hashlib
import json
from pathlib import Path

from audit_target_collisions import audit as census
from audit_collision_engine_sources import audit as engine_audit
from build_collision_evidence import build as evidence

ROOT=Path(__file__).resolve().parents[1]


def replay(cache, output, download=False):
    pinned=json.loads((ROOT/'results/collision_audit/input_manifest.json').read_text())
    if download:
        import public_collision_fetch as fetcher
        fetcher.CACHE=cache
    for record in pinned['requests']:
        name=record['filename']
        if Path(name).name!=name:
            raise ValueError('Manifest must contain plain filenames')
        path=cache/name
        if not path.exists():
            if not download:
                raise FileNotFoundError(f'{name}: supply preserved public cache, or explicitly use --download')
            fetcher.fetch(record['url'],name,body=record.get('request_body'),cap=record['received_bytes'])
        if hashlib.sha256(path.read_bytes()).hexdigest()!=record['sha256']:
            raise ValueError(f'Pinned snapshot changed: {name}. Do not claim this as a reproduction.')
    # Use original snapshot metadata for byte-identical analysis outputs. Leave
    # actual acquisition receipts and their cumulative byte budget untouched.
    census(cache,output,['2024-06','2024-12','2026-06','2026-10'],manifest=pinned)
    engine_audit(cache,output,manifest=pinned)
    evidence(cache,output,manifest=pinned)
    expected=json.loads((ROOT/'results/collision_audit/expected_outputs.json').read_text())
    observed={name:hashlib.sha256((output/name).read_bytes()).hexdigest() for name in expected}
    if observed!=expected:
        mismatched=[name for name in expected if observed[name]!=expected[name]]
        raise ValueError(f'Results differ from recorded audit: {mismatched}')
    print(json.dumps({'reproduced':True,'verified_output_files':len(expected),
                      'verified_public_inputs':len(pinned['requests'])}))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache',type=Path,default=ROOT/'work/public-collision-audit')
    parser.add_argument('--output',type=Path,default=ROOT/'work/collision-replay-results')
    parser.add_argument('--download',action='store_true')
    args=parser.parse_args()
    replay(args.cache,args.output,args.download)
