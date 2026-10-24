"""Independent offline census; deliberately does not import the production audit."""
from collections import defaultdict, Counter
import csv
import gzip
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import argparse
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--cache', type=Path, default=Path(__file__).resolve().parents[1] / 'work/public-collision-audit')
BASE = parser.parse_args().cache

def tab(name):
    with (BASE / name).open(newline='', encoding='utf-8-sig') as stream:
        return list(csv.DictReader((s for s in stream if not s.startswith('#') and not s.startswith('"#')), delimiter='\t'))

def true(s):
    return (s or '').casefold().strip() in ('yes', 'true', 't', '1')

expected = defaultdict(set)
expected_symbols = defaultdict(set)
expected_entrez = defaultdict(set)
for record in tab('gtop_targets.tsv'):
    identifier = record['Target id']
    expected[identifier].update('hgnc:' + v.replace('HGNC:', '') for v in record['HGNC id'].split('|') if v)
    expected_symbols[identifier].update(record['HGNC symbol'].split('|'))
    expected_entrez[identifier].update(record['Human Entrez Gene'].split('|'))

ncbi = {}
ids_for_symbol = {}
with gzip.open(BASE / 'human_gene_info.gz', 'rt') as stream:
    for record in csv.DictReader(stream, delimiter='\t'):
        if record['#tax_id'] != '9606':
            continue
        ncbi[record['GeneID']] = {'hgnc:' + x.rsplit(':', 1)[1] for x in record['dbXrefs'].split('|') if x.startswith('HGNC:')}
        ids_for_symbol[record['Symbol']] = record['GeneID']

source = set()
for record in tab('gtop_interactions.tsv'):
    if record['Target Species'] == 'Human' and not record['Target Ligand']:
        source.add((record['Target ID'], record['Ligand ID'], record['Type'].casefold()))
ligands = {r['Ligand ID']:r for r in tab('gtop_ligands.tsv')}
results = {}
by_gene = defaultdict(list)
all_current = []
for release in ('2024-12','2026-06','2026-10'):
    all_rows = tab(f'dgidb_{release}_interactions.tsv')
    rows = [r for r in all_rows if r['interaction_source_db_name'] == 'GuideToPharmacology']
    stats = Counter()
    proven = []
    unknown = []
    for record in rows:
        prefix, number = record['gene_claim_name'].split(':',1)
        observed = record['gene_concept_id']
        if prefix not in ('NCBIGENE','IUPHAR.RECEPTOR'):
            stats['unsupported_namespace'] += 1
        elif not expected[number]:
            stats['source_unresolved'] += 1
        elif not observed:
            stats['downstream_unresolved'] += 1
        elif observed in expected[number]:
            stats['agrees'] += 1
        elif prefix == 'NCBIGENE' and observed in ncbi.get(number, set()):
            stats['exact_numeric_collision'] += 1
            proven.append(record)
            if release == '2024-12':
                by_gene[observed].append(record)
        else:
            stats['other_conflict'] += 1
            unknown.append(record)
    results[release] = {
        'total_rows':len(all_rows),'gtop_rows':len(rows),'stats':dict(stats),
        'proven_rows':len(proven),
        'proven_source_ligand_pairs':len({(r['gene_concept_id'],r['drug_claim_name']) for r in proven}),
        'proven_normalized_drug_pairs':len({(r['gene_concept_id'],r['drug_concept_id']) for r in proven if r['drug_concept_id']}),
        'proven_genes':len({r['gene_concept_id'] for r in proven}),
        'proven_current_source_match':sum((r['gene_claim_name'].split(':')[1],r['drug_claim_name'].split(':')[1],r['interaction_types'].casefold()) in source for r in proven),
        'proven_dgidb_approved_flag_rows':sum(true(r['drug_is_approved']) for r in proven),
        'unknown_genes':dict(Counter(r['gene_name'] for r in unknown)),
    }
    if release == '2026-10':
        all_current = [{
            'actual_gene':r['gene_name'], 'actual_hgnc':r['gene_concept_id'],
            'actual_gene_ncbi':ids_for_symbol.get(r['gene_name']),
            'source_target_id':r['gene_claim_name'].split(':')[1],
            'ncbi_lookalike_hgnc':sorted(ncbi.get(r['gene_claim_name'].split(':')[1],set())),
            'intended_gene':sorted(expected_symbols[r['gene_claim_name'].split(':')[1]]),
            'intended_hgnc':sorted(expected[r['gene_claim_name'].split(':')[1]]),
            'drug':r['drug_name'], 'ligand':r['drug_claim_name'],
        } for r in unknown]

root = ET.parse(BASE/'orphanet_genes.xml').getroot()
diseases=[]
for disease in root.findall('.//Disorder'):
    affected = set()
    for association in disease.findall('DisorderGeneAssociationList/DisorderGeneAssociation'):
        if association.findtext('DisorderGeneAssociationType/Name') != 'Disease-causing germline mutation(s) in':
            continue
        if association.findtext('DisorderGeneAssociationStatus/Name') != 'Assessed':
            continue
        refs = association.findall('Gene/ExternalReferenceList/ExternalReference')
        ids = {'hgnc:' + r.findtext('Reference').replace('HGNC:', '') for r in refs if r.findtext('Source') == 'HGNC'}
        affected.update(ids & set(by_gene))
    bad=[row for gene in affected for row in by_gene[gene]]
    if bad:
        approved={row['drug_claim_name'].split(':')[1] for row in bad
                  if true(ligands.get(row['drug_claim_name'].split(':')[1],{}).get('Approved'))
                  and not true(ligands.get(row['drug_claim_name'].split(':')[1],{}).get('Withdrawn'))}
        diseases.append({'orpha':disease.findtext('OrphaCode'),'name':disease.findtext('Name'),
                         'group':disease.findtext('DisorderGroup/Name'),
                         'gene_ligand_pairs':len({(r['gene_concept_id'],r['drug_claim_name']) for r in bad}),
                         'approved_ligands':len(approved),
                         'genes':sorted({r['gene_name'] for r in bad})})
diseases.sort(key=lambda d:(-d['approved_ligands'],-d['gene_ligand_pairs'],d['orpha']))
results['orphanet']={'date':root.attrib['date'],'entities':len(diseases),
                      'with_approved_ligands':sum(x['approved_ligands']>0 for x in diseases),
                      'groups':dict(Counter(x['group'] for x in diseases)),
                      'top_15':diseases[:15]}
results['current_conflicts']=all_current
results['target_hgnc_cardinalities']=dict(Counter(len(v) for v in expected.values()))
results['input_hashes']={name:hashlib.sha256((BASE/name).read_bytes()).hexdigest() for name in
    ['gtop_targets.tsv','gtop_interactions.tsv','gtop_ligands.tsv','human_gene_info.gz',
     'dgidb_2024-12_interactions.tsv','dgidb_2026-06_interactions.tsv','dgidb_2026-10_interactions.tsv','orphanet_genes.xml']}
print(json.dumps(results,indent=2))
