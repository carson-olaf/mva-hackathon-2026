"""Audit current engine sources against the observed historical collision set.

This is not a claim that all assertions in these databases are correct.
"""
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re

from audit_target_collisions import hgnc, target_map, tsv


def audit(cache, output, manifest=None):
    if manifest is None:
        manifest=json.loads((cache/'manifest.json').read_text())
    records={r['filename']:r for r in manifest['requests'] if r['complete']}
    used=[]
    def verify(name):
        raw=(cache/name).read_bytes()
        if hashlib.sha256(raw).hexdigest()!=records[name]['sha256']:
            raise ValueError(f'Hash mismatch: {name}')
        used.append(records[name])
        return raw
    def read(name):
        return json.loads(verify(name))
    for name in ['gtop_ligands.tsv','gtop_targets.tsv','human_gene_info.gz']:
        verify(name)
    mechanism_files=sorted(cache.glob('chembl_mechanisms_*.json'))
    pages=[read(f.name) for f in mechanism_files]
    mechanisms=[r for page in pages for r in page['mechanisms']]
    expected={page['page_meta']['total_count'] for page in pages}
    if expected!={len(mechanisms)} or len({r['mec_id'] for r in mechanisms})!=len(mechanisms):
        raise ValueError('Incomplete, overlapping or changing ChEMBL mechanism collection')
    targets=[r for f in sorted(cache.glob('chembl_targets_*.json')) for r in read(f.name)['targets']]
    target_by_id={r['target_chembl_id']:r for r in targets}
    wanted={r['target_chembl_id'] for r in mechanisms if r['target_chembl_id']}
    if set(target_by_id)!=wanted or len(target_by_id)!=len(targets):
        raise ValueError('ChEMBL target set incomplete or duplicated')
    with (output/'identity_conflicts.csv').open() as handle:
        bad=[r for r in csv.DictReader(handle) if r['snapshot']=='2024-12' and r['classification']=='confirmed_numeric_namespace_collision']
    ligands={r['Ligand ID']:r for r in tsv(cache/'gtop_ligands.tsv')}
    comparable=defaultdict(list)
    corrected_comparable=set()
    mapped_source_pairs=set()
    for row in bad:
        for chembl in ligands.get(row['ligand_id'],{}).get('ChEMBL ID','').split('|'):
            if chembl:
                comparable[(row['wrong_hgnc'],chembl)].append(row)
                corrected_comparable.update((g,chembl) for g in row['correct_hgnc'].split('|') if g)
                mapped_source_pairs.add((row['wrong_hgnc'],row['ligand_id']))
    overlaps=set()
    corrected_overlaps=set()
    human=0
    human_phase4=0
    hgnc_resolved=0
    source_map=target_map(tsv(cache/'gtop_targets.tsv'))
    by_uniprot=defaultdict(set)
    for target in source_map.values():
        for protein in target['uniprot']:
            by_uniprot[protein].update(target['hgnc'])
    protein_checks=set()
    protein_conflicts=set()
    for target in targets:
        if target['tax_id']!=9606:
            continue
        for component in target['target_components']:
            gids={hgnc(x['xref_id']) for x in component['target_component_xrefs'] if x['xref_src_db']=='HGNC'}
            accession=component['accession']
            if gids and by_uniprot.get(accession):
                protein_checks.add(accession)
                if not gids.intersection(by_uniprot[accession]):
                    protein_conflicts.add(accession)
    for row in mechanisms:
        target=target_by_id.get(row['target_chembl_id'],{})
        if target.get('tax_id')!=9606:
            continue
        human+=1
        human_phase4+=row['max_phase']==4
        gids={hgnc(x['xref_id']) for c in target['target_components'] for x in c['target_component_xrefs'] if x['xref_src_db']=='HGNC'}
        hgnc_resolved+=bool(gids)
        for gene in gids:
            for compound in {row['molecule_chembl_id'],row['parent_molecule_chembl_id']}:
                if (gene,compound) in comparable:
                    overlaps.add((gene,compound,row['mec_id']))
                if (gene,compound) in corrected_comparable:
                    corrected_overlaps.add((gene,compound))
    # Check every requested Ensembl ID against the actual returned ID. Symbols
    # are a supplementary comparison; a renamed symbol alone is not a collision.
    affected_entrez={r['claim_id'].partition(':')[2] for r in bad}
    request_map=[]
    with gzip.open(cache/'human_gene_info.gz','rt') as handle:
        for gene in csv.DictReader(handle,delimiter='\t'):
            if gene['#tax_id']!='9606' or gene['GeneID'] not in affected_entrez:
                continue
            ensembl={x.partition(':')[2] for x in gene['dbXrefs'].split('|') if x.startswith('Ensembl:')}
            if len(ensembl)==1:
                request_map.append([gene['GeneID'],gene['Symbol'],next(iter(ensembl))])
    request_map.sort(key=lambda r:int(r[0]))
    ot_counts=Counter()
    ot_disagreements=[]
    ot_protein_checks=set()
    ot_protein_disagreements=[]
    chembl_protein_symbols=defaultdict(set)
    for target in targets:
        if target['tax_id']==9606:
            for c in target['target_components']:
                chembl_protein_symbols[c['accession']].update(s['component_synonym'] for s in c['target_component_synonyms'] if s['syn_type']=='GENE_SYMBOL')
    for offset in range(0,len(request_map),100):
        name=f'opentargets_targets_{offset}.json'
        response=read(name)
        if response.get('errors'):
            raise ValueError(f'Open Targets errors: {name}')
        for index,(entrez,symbol,ensembl) in enumerate(request_map[offset:offset+100]):
            # Also prove that the stored request used the expected identifier.
            if f'g{index}:target(ensemblId:"{ensembl}")' not in records[name]['request_body']['query']:
                raise ValueError('Open Targets request map mismatch')
            ot_counts['requested']+=1
            row=response['data'][f'g{index}']
            if row is None:
                ot_counts['unresolved']+=1
                continue
            ot_counts['resolved']+=1
            if row['id']!=ensembl:
                ot_disagreements.append({'kind':'ensembl_id','requested':ensembl,'returned':row['id']})
            if row['approvedSymbol']!=symbol:
                ot_disagreements.append({'kind':'symbol','requested':symbol,'returned':row['approvedSymbol'],'ensembl':ensembl})
            for protein in row['proteinIds']:
                if protein['source']=='uniprot_swissprot' and protein['id'] in chembl_protein_symbols:
                    ot_protein_checks.add((ensembl,protein['id']))
                    if row['approvedSymbol'] not in chembl_protein_symbols[protein['id']]:
                        ot_protein_disagreements.append({'ensembl':ensembl,'uniprot':protein['id'],'ot_symbol':row['approvedSymbol'],'chembl_symbols':sorted(chembl_protein_symbols[protein['id']])})
    summary={
        'chembl_status':read('chembl_status.json'),
        'chembl_mechanism_records':len(mechanisms),
        'chembl_unique_targets':len(targets),
        'human_mechanism_records':human,
        'human_phase4_mechanism_records':human_phase4,
        'human_mechanism_records_with_hgnc_target':hgnc_resolved,
        'historical_wrong_gene_ligand_pairs_with_chembl_ligand_map':len(mapped_source_pairs),
        'historical_wrong_gene_chembl_compound_pairs_after_mapping':len(comparable),
        'historical_wrong_pair_overlaps_in_current_mechanisms':sorted(overlaps),
        'corrected_source_counterpart_pairs_in_current_mechanisms':len(corrected_overlaps),
        'chembl_gtop_shared_protein_identity_checks':len(protein_checks),
        'chembl_gtop_protein_identity_conflicts':sorted(protein_conflicts),
        'opentargets':dict(ot_counts),
        'opentargets_identifier_or_symbol_disagreements':ot_disagreements,
        'opentargets_chembl_protein_checks':len(ot_protein_checks),
        'opentargets_chembl_symbol_disagreements':ot_protein_disagreements,
        'reactome_scope':'Engine uses Reactome UniProt pathway membership, not a direct drug-gene assertion source. No exhaustive Reactome content audit claimed.',
        'limits':['No ChEMBL match cannot exclude unannotated biological effects.',
                  'Only historical pairs with a GtoP ligand-to-ChEMBL cross-reference can be searched.',
                  'Other database tables and error classes are outside this census.',
                  'Symbol disagreements need alias review; they are not automatically numeric collisions.'],
        'inputs':used,
    }
    (output/'engine_source_audit.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k!='inputs'},indent=2))


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache',type=Path,default=Path('work/public-collision-audit'))
    parser.add_argument('--output',type=Path,default=Path('results/collision_audit'))
    args=parser.parse_args()
    audit(args.cache,args.output)
