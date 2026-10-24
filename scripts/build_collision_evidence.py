"""Generate named-case and live-sentinel evidence from verified public inputs."""
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path

from audit_target_collisions import truth, write_csv

# Selection is editorial; all counts and identity evidence are computed.
# PMID 16244177 was traced separately because its GtoP interaction citation
# field is blank. The other pharmacology PMIDs occur in the source records.
CASES = [
    ('GM2A', '7335', '309246', '32106752'),
    ('GATM', '6780', '3337', '15102696'),
    ('ABCD1', '958', '139396', '8301582'),
    ('AARS1', '320', '228174', '27085897'),
    ('ABCA4', '485', '827', '34355529'),
    ('CPT2', '1909', '228302', '29031075'),
    ('BUB1B', '2439', '1052', '16244177'),
]


def build(cache, output, manifest=None):
    if manifest is None:
        manifest=json.loads((cache/'manifest.json').read_text())
    records={r['filename']:r for r in manifest['requests'] if r['complete']}
    used=[]
    def read(name):
        raw=(cache/name).read_bytes()
        if hashlib.sha256(raw).hexdigest()!=records[name]['sha256']:
            raise ValueError(f'Hash mismatch: {name}')
        used.append(records[name])
        return json.loads(raw)
    publications={r['id']:r for r in read('case_publications.json')['resultList']['result'] if r['source']=='MED'}
    with (output/'identity_conflicts.csv').open() as handle:
        historic=[r for r in csv.DictReader(handle) if r['snapshot']=='2024-12' and r['classification']=='confirmed_numeric_namespace_collision']
    with (output/'rare_disease_exposure.csv').open() as handle:
        diseases={r['orpha']:r for r in csv.DictReader(handle)}
    cases=[]
    for gene,ligand,orpha,pmid in CASES:
        candidates=[r for r in historic if r['wrong_gene']==gene and r['ligand_id']==ligand]
        if len(candidates)!=1:
            raise ValueError(f'Expected one observed case: {gene}/{ligand}')
        case=candidates[0].copy()
        if gene!='BUB1B' and pmid not in case['source_pmids'].split('|'):
            raise ValueError(f'Case PMID is not in the source evidence: {gene}')
        pub=publications[pmid]
        case.update(
            orpha=orpha, disease=diseases[orpha]['disease'],
            wrong_source_pairs_for_gene=len({r['ligand_id'] for r in historic if r['wrong_gene']==gene}),
            disease_approved_flag_normalized_drugs=diseases[orpha]['unique_approved_flag_normalized_drugs'],
            opened_pharmacology_pmid=pmid, opened_pharmacology_doi=pub.get('doi',''),
            opened_pharmacology_title=pub['title'],
            opened_scope='Abstract; PMID32106752 also PubMed Figure1 legend',
        )
        cases.append(case)
    write_csv(output/'named_cases.csv',cases)
    live=read('dgidb_live_sentinels.json')
    if live.get('errors'):
        raise ValueError('DGIdb sentinel query returned errors')
    query=records['dgidb_live_sentinels.json']['request_body']['query']
    checked=[]
    for alias,gid in [('g0','1149'),('g1','1148'),('g2','32041'),('g3','21160'),('g4','6943'),('g5','8772')]:
        if f'{alias}:gene(conceptId:"hgnc:{gid}")' not in query:
            raise ValueError('Sentinel request identity mismatch')
        row=live['data'][alias]
        gtop=[] if row is None else [i for i in row['interactions'] if any(s['fullName']=='Guide to Pharmacology' for s in i['sources'])]
        checked.append({'hgnc':f'hgnc:{gid}','name':None if row is None else row['name'],
                        'total_interactions':None if row is None else len(row['interactions']),
                        'gtop_interactions':gtop})
    with (output/'current_identity_conflicts.csv').open() as handle:
        conflicts=list(csv.DictReader(handle))
    confirmed=[]
    for row in conflicts:
        sentinel=next(s for s in checked if s['hgnc']==row['wrong_hgnc'])
        matches=[i for i in sentinel['gtop_interactions'] if i['drug']['conceptId']==row['drug_concept_id']]
        if len(matches)!=1:
            raise ValueError(f'Live/export disagreement: {row["wrong_gene"]}/{row["drug"]}')
        confirmed.append({'assigned_gene':row['wrong_gene'],'correct_gene':row['correct_gene'],
                          'source_target_id':row['gtop_target_id'],'drug':row['drug'],
                          'drug_concept_id':row['drug_concept_id']})
    result={'retrieved_utc':records['dgidb_live_sentinels.json']['retrieved_utc'],
            'confirmed_live_identity_conflicts':confirmed,'sentinels':checked,
            'scope':'Only listed Guide to Pharmacology edges are adjudicated; other returned source edges are not counted as errors.',
            'inputs':used}
    (output/'live_api_evidence.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'named_cases':len(cases),'live_conflicts':len(confirmed)}))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache',type=Path,default=Path('work/public-collision-audit'))
    parser.add_argument('--output',type=Path,default=Path('results/collision_audit'))
    args=parser.parse_args()
    build(args.cache,args.output)
