"""Offline public-data census of GtoP target IDs mislabelled as NCBI Gene IDs.

This reads only the explicitly named public cache. It never reads patient data.
Counts refer to observed export rows, not all possible numeric ID collisions.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET


def tsv(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader((line for line in handle if not line.lstrip('"').startswith('#')), delimiter='\t'))


def dgidb_rows(path):
    rows = tsv(path)
    if rows and not {'gene_concept_id','interaction_source_db_name','interaction_source_db_version'}.issubset(rows[0]):
        raise ValueError('Expected DGIdb v5 export schema; legacy Entrez claims have different semantics')
    for row in rows:
        if 'interaction_types' not in row:
            row['interaction_types'] = row.pop('interaction_type')
        if 'drug_is_approved' not in row:
            row['drug_is_approved'] = row.pop('approved').lower()
        if row['interaction_types'] == 'NULL':
            row['interaction_types'] = ''
    return rows


def hgnc(value):
    value = value.removeprefix('HGNC:').removeprefix('hgnc:')
    if not re.fullmatch(r'\d+', value):
        raise ValueError(f'Not a human HGNC identifier: {value}')
    return 'hgnc:' + value


def truth(value):
    return str(value).strip().lower() in {'t','true','yes','1'}


def target_map(rows):
    targets = {}
    for row in rows:
        key = row['Target id']
        target = targets.setdefault(key, {'hgnc': set(), 'entrez': set(), 'symbols': set(), 'uniprot': set()})
        for field, name in [('HGNC id','hgnc'), ('Human Entrez Gene','entrez'), ('HGNC symbol','symbols'), ('Human SwissProt','uniprot')]:
            for value in row[field].split('|'):
                if not value:
                    continue
                if name == 'hgnc':
                    try:
                        target[name].add(hgnc(value))
                    except ValueError:
                        # Preserve the anomaly for reporting; never invent an
                        # HGNC ID by prefixing an MGI or other foreign ID.
                        target.setdefault('invalid_hgnc_values',set()).add(value)
                else:
                    target[name].add(value)
    return targets


def classify(row, targets, ncbi):
    """Identity disagreement is not automatically evidence of the numeric bug."""
    claim = row['gene_claim_name']
    namespace, _, number = claim.partition(':')
    if namespace not in {'NCBIGENE', 'IUPHAR.RECEPTOR'}:
        return 'unsupported_claim_namespace'
    target = targets.get(number)
    if not target or not target['hgnc']:
        return 'unresolved_source_target'
    if not row['gene_concept_id']:
        return 'unresolved_downstream_gene'
    if row['gene_concept_id'] in target['hgnc']:
        return 'identity_agrees'
    if namespace == 'NCBIGENE' and row['gene_concept_id'] in ncbi.get(number, set()):
        return 'confirmed_numeric_namespace_collision'
    return 'identity_conflict_other_or_unresolved_cause'


def write_csv(path, rows, fields=None):
    if not rows and not fields:
        # An empty result must replace a previous run, not leave stale findings.
        path.write_text('')
        return
    fields = fields or list(rows[0])
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def count_edges(rows):
    return {
        'rows': len(rows),
        'distinct_gene_ligand_pairs': len({(r['gene_concept_id'], r['drug_claim_name']) for r in rows}),
        'distinct_gene_normalized_drug_pairs': len({(r['gene_concept_id'], r['drug_concept_id']) for r in rows if r['drug_concept_id']}),
        'distinct_gene_ligand_type_triples': len({(r['gene_concept_id'], r['drug_claim_name'], r['interaction_types']) for r in rows}),
        'distinct_genes': len({r['gene_concept_id'] for r in rows}),
        'distinct_ligands': len({r['drug_claim_name'] for r in rows}),
    }


def audit(cache, output, tags, exposure_snapshot='2024-12', manifest=None):
    output.mkdir(parents=True, exist_ok=True)
    if manifest is None:
        manifest = json.loads((cache/'manifest.json').read_text())
    required = ['gtop_targets.tsv','gtop_interactions.tsv','gtop_ligands.tsv','human_gene_info.gz','orphanet_genes.xml']
    required += [f'dgidb_{tag}_interactions.tsv' for tag in tags]
    inputs = []
    for filename in required:
        records = [r for r in manifest['requests'] if r['filename'] == filename and r['complete']]
        if not records:
            raise ValueError(f'Missing completed public input: {filename}')
        record = records[-1]
        if hashlib.sha256((cache/filename).read_bytes()).hexdigest() != record['sha256']:
            raise ValueError(f'Hash mismatch: {filename}')
        inputs.append(record)
    targets = target_map(tsv(cache/'gtop_targets.tsv'))
    ncbi = {}
    with gzip.open(cache/'human_gene_info.gz', 'rt') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            if row['#tax_id'] == '9606':
                ncbi[row['GeneID']] = {hgnc(x.split(':')[-1]) for x in row['dbXrefs'].split('|') if x.startswith('HGNC:')}
    source_edges = defaultdict(list)
    source_pairs = defaultdict(list)
    for row in tsv(cache/'gtop_interactions.tsv'):
        if row['Target Species'] == 'Human' and not row['Target Ligand']:
            source_edges[(row['Target ID'], row['Ligand ID'], row['Type'].lower())].append(row)
            source_pairs[(row['Target ID'], row['Ligand ID'])].append(row)
    ligands = {r['Ligand ID']: r for r in tsv(cache/'gtop_ligands.tsv')}
    snapshots = {}
    conflicts = []
    observed = {}
    for tag in tags:
        all_rows = dgidb_rows(cache/f'dgidb_{tag}_interactions.tsv')
        rows = [r for r in all_rows if r['interaction_source_db_name'] == 'GuideToPharmacology']
        observed[tag] = rows
        statuses = Counter()
        confirmed = []
        backed = []
        pair_backed = []
        for row in rows:
            status = classify(row, targets, ncbi)
            statuses[status] += 1
            if status == 'identity_agrees':
                continue
            tid = row['gene_claim_name'].partition(':')[2]
            lid = row['drug_claim_name'].partition(':')[2]
            target = targets.get(tid, {'hgnc':set(),'symbols':set(),'entrez':set(),'uniprot':set()})
            evidence = source_edges.get((tid,lid,row['interaction_types'].lower()), [])
            pair_evidence = source_pairs.get((tid,lid), [])
            ligand = ligands.get(lid, {})
            record = {
                'snapshot':tag, 'classification':status,
                'wrong_gene':row['gene_name'], 'wrong_hgnc':row['gene_concept_id'],
                'claim_id':row['gene_claim_name'], 'gtop_target_id':tid,
                'correct_gene':'|'.join(sorted(target['symbols'])),
                'correct_hgnc':'|'.join(sorted(target['hgnc'])),
                'correct_entrez':'|'.join(sorted(target['entrez'])),
                'correct_uniprot':'|'.join(sorted(target['uniprot'])),
                'ligand_id':lid, 'drug':row['drug_name'],
                'source_ligand_name':ligand.get('Name',''),
                'drug_concept_id':row['drug_concept_id'], 'interaction_type':row['interaction_types'],
                'dgidb_approved_flag':row['drug_is_approved'],
                'gtop_approved_flag':ligand.get('Approved',''),
                'gtop_withdrawn_flag':ligand.get('Withdrawn',''),
                'dgidb_source_version':row['interaction_source_db_version'],
                'current_source_same_target_ligand_type':bool(evidence),
                'current_source_same_target_ligand':bool(pair_evidence),
                'current_source_types':'|'.join(sorted({e['Type'] for e in pair_evidence})),
                'source_pmids':'|'.join(sorted({p for e in pair_evidence for p in e['PubMed ID'].split('|') if p})),
                'source_target_url':f'https://www.guidetopharmacology.org/GRAC/ObjectDisplayForward?objectId={tid}',
            }
            conflicts.append(record)
            if status == 'confirmed_numeric_namespace_collision':
                confirmed.append(row)
                if evidence:
                    backed.append(row)
                if pair_evidence:
                    pair_backed.append(row)
        snapshots[tag] = {
            'all_source_rows':len(all_rows), 'gtop_source_rows':len(rows),
            'source_versions':dict(Counter(r['interaction_source_db_version'] for r in rows)),
            'classifications':dict(statuses),
            'confirmed_numeric_collisions':count_edges(confirmed),
            'confirmed_with_current_source_edge':count_edges(backed),
            'confirmed_with_current_source_pair':count_edges(pair_backed),
            'confirmed_dgidb_approved_flag':count_edges([r for r in confirmed if r['drug_is_approved']=='true']),
        }
    # A release comparison supplies independent observation of repaired assignments.
    observed_index = {tag: {
        (r['gene_concept_id'],r['gene_claim_name'].partition(':')[2],
         r['drug_claim_name'].partition(':')[2],r['interaction_types']) for r in rows
    } for tag, rows in observed.items()}
    for record in conflicts:
        expected = set(record['correct_hgnc'].split('|'))
        for tag in [t for t in tags if t >= '2026-06']:
            record[f'correct_pair_present_{tag}'] = any(
                (gene,record['gtop_target_id'],record['ligand_id'],record['interaction_type']) in observed_index[tag]
                for gene in expected
            )
    write_csv(output/'identity_conflicts.csv', conflicts)
    current = [r for r in conflicts if r['snapshot'] == tags[-1]]
    write_csv(output/'current_identity_conflicts.csv', current, list(conflicts[0]) if conflicts else None)
    # Rank potential misleading leads, not clinical harm: assessed germline
    # disease links and the number of distinct misplaced approved-flag ligands.
    historic = [r for r in conflicts if r['snapshot']==exposure_snapshot and r['classification']=='confirmed_numeric_namespace_collision']
    by_gene = defaultdict(list)
    for row in historic:
        by_gene[row['wrong_hgnc']].append(row)
    root = ET.parse(cache/'orphanet_genes.xml').getroot()
    diseases = []
    associations = []
    for disorder in root.iter('Disorder'):
        ids = set()
        for assoc in disorder.findall('./DisorderGeneAssociationList/DisorderGeneAssociation'):
            kind = assoc.findtext('./DisorderGeneAssociationType/Name','')
            status = assoc.findtext('./DisorderGeneAssociationStatus/Name','')
            if kind != 'Disease-causing germline mutation(s) in' or status != 'Assessed':
                continue
            for ref in assoc.findall('./Gene/ExternalReferenceList/ExternalReference'):
                if ref.findtext('Source')=='HGNC':
                    gid=hgnc(ref.findtext('Reference'))
                    ids.add(gid)
                    if gid in by_gene:
                        associations.append({'orpha':disorder.findtext('OrphaCode'), 'disease':disorder.findtext('Name'),
                            'hgnc':gid, 'gene':assoc.findtext('./Gene/Symbol'),
                            'validation_sources':assoc.findtext('SourceOfValidation','')})
        hits = [r for gid in ids for r in by_gene.get(gid, [])]
        if hits:
            approved = [r for r in hits if truth(r['gtop_approved_flag']) and not truth(r['gtop_withdrawn_flag'])]
            diseases.append({'orpha':disorder.findtext('OrphaCode'), 'disease':disorder.findtext('Name'),
                'disorder_group':disorder.findtext('./DisorderGroup/Name'),
                'affected_genes':'|'.join(sorted({r['wrong_gene'] for r in hits})),
                'unique_wrong_gene_ligand_pairs':len({(r['wrong_hgnc'],r['ligand_id']) for r in hits}),
                'unique_approved_flag_gene_ligand_pairs':len({(r['wrong_hgnc'],r['ligand_id']) for r in approved}),
                'unique_approved_flag_ligands':len({r['ligand_id'] for r in approved}),
                'unique_approved_flag_normalized_drugs':len({r['drug_concept_id'] for r in approved if r['drug_concept_id']}),
                'approved_flag_pairs_without_normalized_drug':sum(not r['drug_concept_id'] for r in approved),
                'approved_flag_ligand_names':'|'.join(sorted({r['drug'] for r in approved})),
            })
    diseases.sort(key=lambda r:(-r['unique_approved_flag_normalized_drugs'],-r['unique_wrong_gene_ligand_pairs'],r['orpha']))
    write_csv(output/'rare_disease_exposure.csv',diseases)
    write_csv(output/'rare_disease_gene_evidence.csv',associations)
    summary = {
        'scope':'Observed DGIdb GuideToPharmacology claims, independently reconciled with GtoP and human NCBI-HGNC IDs',
        'snapshots':snapshots,
        'invalid_hgnc_values_in_gtop_targets': {k:sorted(v['invalid_hgnc_values']) for k,v in targets.items() if v.get('invalid_hgnc_values')},
        'rare_disease_map_date':root.attrib.get('date'),
        'rare_disease_exposure_snapshot':exposure_snapshot,
        'rare_disease_association_filter':'Assessed; Disease-causing germline mutation(s) in; exact HGNC ID',
        'exposed_orpha_entities':len(diseases),
        'exposed_orpha_entities_with_approved_flag_ligand':sum(r['unique_approved_flag_ligands']>0 for r in diseases),
        'exposed_orpha_entities_by_group':dict(Counter(r['disorder_group'] for r in diseases)),
        'exposed_orpha_genes':len({r['hgnc'] for r in associations}),
        'disease_ranking':'Distinct normalized drug concepts with current GtoP approved/nonwithdrawn flags; source-ligand counts retained separately',
        'limits':['Current GtoP evidence is not an archived 2024 source file.',
                  'Approval flags are database annotations, not independent current regulator checks.',
                  'Other identity conflicts are not counted as confirmed numeric namespace collisions.',
                  'Disease exposure describes erroneous search leads, not measured patient harm.',
                  'These counts do not constitute a census of every database or all error classes.'],
        'inputs':inputs,
    }
    (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    return summary


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache',type=Path,default=Path('work/public-collision-audit'))
    parser.add_argument('--output',type=Path,default=Path('results/collision_audit'))
    parser.add_argument('--snapshots',nargs='+',default=['2024-06','2024-12','2026-06','2026-10'])
    parser.add_argument('--exposure-snapshot',default='2024-12')
    args=parser.parse_args()
    result=audit(args.cache,args.output,args.snapshots,args.exposure_snapshot)
    print(json.dumps({k:v for k,v in result.items() if k!='inputs'},indent=2))
