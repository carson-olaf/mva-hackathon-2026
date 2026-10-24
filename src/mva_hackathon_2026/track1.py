"""Phenotype-first, gene-blind discovery. Patient-scale records stay in work/."""
import argparse
import csv
import hashlib
import itertools
import json
import math
import subprocess
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import quote
from .public_data import PublicCache

HPO = ['HP:0002859','HP:0000121','HP:0004322','HP:0001508','HP:0003202','HP:0001622','HP:0001518']
# The submitted run also used one family-history HPO term, withheld from the public copy to protect the family's privacy.
HPO_URL = 'https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/phenotype_to_genes.txt'
BASE = 'https://rest.ensembl.org'
CODE = dict(zip((''.join(c) for c in itertools.product('TCAG', repeat=3)), 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))
IMPACT = {'stop_gained':6,'frameshift_variant':6,'inframe_indel':3,'missense_variant':2,'stop_lost':3,'synonymous_variant':0}


def revcomp(seq):
    return seq.translate(str.maketrans('ACGTN','TGCAN'))[::-1]


def translate(seq):
    return ''.join(CODE.get(seq[i:i+3],'X') for i in range(0,len(seq)-2,3))


def phenotype_prior(text, terms):
    by_term, genes = defaultdict(set), set()
    for line in text.splitlines():
        fields=line.split('\t')
        if len(fields)>=4 and fields[0].startswith('HP:'):
            by_term[fields[0]].add(fields[3]); genes.add(fields[3])
    weights={t:math.log((len(genes)+1)/(len(by_term[t])+1)) for t in terms}
    prior=defaultdict(float)
    for term in terms:
        for gene in by_term[term]: prior[gene]+=weights[term]
    return dict(prior),{t:{'gene_count':len(by_term[t]),'weight':weights[t]} for t in terms},by_term


def gene_model(cache,symbol):
    gene=cache.get(f'{BASE}/lookup/symbol/homo_sapiens/{symbol}',params={'expand':1})
    transcripts=[t for t in gene.get('Transcript',[]) if t.get('is_canonical') and t.get('Translation')]
    if not transcripts: return None
    tx=transcripts[0]
    seq=cache.get(f"{BASE}/sequence/id/{tx['id']}",params={'type':'cds'})['seq'].upper()
    exons=sorted(tx['Exon'],key=lambda e:e['start'],reverse=tx['strand']==-1)
    start,end=tx['Translation']['start'],tx['Translation']['end']
    positions,exon_numbers=[],[]
    for number,exon in enumerate(exons,1):
        lo,hi=max(start,exon['start']),min(end,exon['end'])
        if lo<=hi:
            segment=list(range(lo,hi+1))
            if tx['strand']==-1: segment.reverse()
            positions.extend(segment);exon_numbers.extend([number]*len(segment))
    if len(seq)!=len(positions): raise ValueError(f'{symbol}: CDS/model length mismatch {len(seq)} vs {len(positions)}')
    return {'symbol':symbol,'chrom':str(gene['seq_region_name']),'start':gene['start'],'end':gene['end'],'transcript':tx['id']+'.'+str(tx['version']),'strand':tx['strand'],'sequence':seq,'positions':positions,'exon_numbers':exon_numbers,'exons':exons}


def annotate(model,pos,ref,alt):
    lookup={p:i for i,p in enumerate(model['positions'])}
    span=[lookup.get(p) for p in range(pos,pos+len(ref))]
    if any(i is None for i in span): return None
    indices=sorted(span)
    if indices!=list(range(indices[0],indices[-1]+1)): return None
    index=indices[0]
    if model['strand']==-1: ref,alt=revcomp(ref),revcomp(alt)
    seq=model['sequence']
    if seq[index:index+len(ref)]!=ref: return {'error':'reference_mismatch'}
    wt,mut=translate(seq),translate(seq[:index]+alt+seq[index+len(ref):])
    aa=index//3
    if (len(alt)-len(ref))%3: effect='frameshift_variant'
    elif len(ref)!=len(alt): effect='inframe_indel'
    elif wt==mut: effect='synonymous_variant'
    elif '*' in mut[:len(wt)-1] and ('*' not in wt[:len(wt)-1] or mut.index('*')<wt.index('*')): effect='stop_gained'
    elif wt.endswith('*') and not mut.endswith('*'): effect='stop_lost'
    else: effect='missense_variant'
    if effect == 'frameshift_variant':
        protein = f'frameshift at/after codon {aa+1}; exact HGVS not computed'
    elif effect == 'inframe_indel':
        protein = f'in-frame indel at/after codon {aa+1}; exact HGVS not computed'
    else:
        protein = f'{wt[aa]}{aa+1}{mut[aa]}'
    return {'consequence':effect,'coding_position':index+1,'protein':protein,'exon':model['exon_numbers'][index],'exon_count':len(model['exons']),'transcript':model['transcript']}


def extract_gene(vcf,model):
    fmt='%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER[\t%GT\t%DP\t%GQ\t%AD]\n'
    result=subprocess.run(['bcftools','query','-r',f"{model['chrom']}:{model['start']}-{model['end']}",'-f',fmt,str(vcf)],check=True,text=True,capture_output=True)
    counts,variants=Counter(),[]
    for line in result.stdout.splitlines():
        chrom,pos,ref,alt,qual,filt,gt,dp,gq,ad=line.split('\t')
        counts['regional_records']+=1
        if filt!='PASS' or ',' in alt or any(s in alt for s in '<>*') or gt not in ('0/1','1/0','0|1','1|0','1/1','1|1'): continue
        if '.' in (qual,dp,gq) or float(qual)<30 or int(dp)<10 or int(gq)<30: continue
        alleles=ad.split(',')
        if len(alleles)!=2 or '.' in alleles: continue
        a,b=map(int,alleles)
        if a+b==0: continue
        het=gt in ('0/1','1/0','0|1','1|0')
        if het and not 0.2<=b/(a+b)<=0.8: continue
        counts['quality_records']+=1
        annotation=annotate(model,int(pos),ref,alt)
        if annotation is None: continue
        if 'error' in annotation: counts['reference_mismatches']+=1;continue
        counts[annotation['consequence']]+=1
        if not IMPACT[annotation['consequence']]: continue
        variants.append({'gene':model['symbol'],'chrom':chrom,'pos':int(pos),'ref':ref,'alt':alt,'gt':gt,'dp':int(dp),'gq':int(gq),'ab':b/(a+b),'het':het,**annotation})
    return variants,counts


def rank_models(variants,prior):
    by_gene=defaultdict(list)
    for v in variants:
        if v.get('max_population_af') is not None and v['max_population_af']>0.01: continue
        by_gene[v['gene']].append(v)
    models=[]
    for gene,calls in by_gene.items():
        for v in calls:
            mode='de_novo_compatible_dominant' if v['het'] else 'recessive_homozygous'
            models.append({'gene':gene,'model':mode,'variants':[v],'score':prior[gene]/3+IMPACT[v['consequence']]+(0 if v['het'] else 3)})
        for pair in itertools.combinations([v for v in calls if v['het']],2):
            models.append({'gene':gene,'model':'recessive_compound_het_unphased','variants':list(pair),'score':prior[gene]/3+sum(IMPACT[v['consequence']] for v in pair)+4})
    return sorted(models,key=lambda x:(-x['score'],x['gene'],x['variants'][0]['pos']))


def vep_shortlist(cache,ranked,maximum=12):
    selected={}
    for model in ranked:
        for v in model['variants']:
            key=(v['chrom'],v['pos'],v['ref'],v['alt'])
            if len(selected)<maximum: selected.setdefault(key,v)
        if len(selected)>=maximum: break
    for v in selected.values():
        if len(v['ref'])!=1 or len(v['alt'])!=1:
            v['population_annotation']='not_queried_indel';continue
        hgvs=f"{v['chrom']}:g.{v['pos']}{v['ref']}>{v['alt']}"
        try:
            data=cache.get(f'{BASE}/vep/human/hgvs/{quote(hgvs,safe="")}',params={'CADD':1,'REVEL':1,'AlphaMissense':1,'canonical':1,'hgvs':1})
            freqs=[]
            by_source={}
            for colocated in data[0].get('colocated_variants',[]):
                for source,f in colocated.get('frequencies',{}).get(v['alt'],{}).items():
                    if isinstance(f,(int,float)):
                        freqs.append(f)
                        by_source[source]=f
            v['max_population_af']=max(freqs) if freqs else None
            v['population_af_by_source']=by_source
            v['population_annotation']='VEP_population_frequencies' if freqs else 'no_frequency_returned_not_proof_of_absence'
            consequence=next((t for t in data[0].get('transcript_consequences',[]) if t.get('transcript_id')==v['transcript'].split('.')[0]),{})
            v['vep']={k:consequence[k] for k in ('hgvsc','hgvsp','consequence_terms','cadd_phred','revel_score','sift_prediction','polyphen_prediction','am_class','am_pathogenicity') if k in consequence}
            if consequence.get('consequence_terms') and v['consequence'] not in consequence['consequence_terms']: v['annotation_discordance']=True
        except Exception as exc:
            v['population_annotation']=f'query_failed_{type(exc).__name__}'
    return len(selected)


def nmd(model,variant):
    exons=model['exons']; cdna_positions=[]
    for exon in exons:
        segment=list(range(exon['start'],exon['end']+1))
        if model['strand']==-1: segment.reverse()
        cdna_positions.extend(segment)
    coding_start=cdna_positions.index(model['positions'][0])
    codon_start=((variant['coding_position']-1)//3)*3
    stop_end_cdna=coding_start+codon_start+3
    last_junction=sum(e['end']-e['start']+1 for e in exons[:-1])
    gap=last_junction-stop_end_cdna
    return {'transcript':model['transcript'],'stop_end_cdna':stop_end_cdna,'last_exon_junction_cdna':last_junction,'bases_after_stop_before_final_junction':gap,'nmd_expected_50_55_rule':gap>55,'caveat':'Prediction only; transcript abundance and alternative splicing require RNA experiments.','reference':'PMID:35052820'}


def run(args):
    work=Path(args.work)
    if not work.resolve().is_relative_to((Path.cwd()/'work').resolve()):
        raise ValueError('Patient intermediates must remain beneath the ignored work/ directory')
    work.mkdir(parents=True,exist_ok=True)
    ignored=subprocess.run(['git','check-ignore','-q',str(work/'track1-candidates.json')])
    if ignored.returncode:
        raise ValueError('Patient intermediate path must be git-ignored')
    output=Path(args.output); output.mkdir(parents=True,exist_ok=True)
    cache=PublicCache(work/'track1-cache')
    terms=args.hpo or HPO
    prior,weights,byterm=phenotype_prior(cache.get(HPO_URL,text=True),terms)
    symbols=sorted(prior,key=lambda g:(-prior[g],g))[:args.gene_limit]
    models,failures={},{}
    def fetch(symbol):
        try: return symbol,gene_model(cache,symbol),None
        except Exception as exc: return symbol,None,str(exc)
    with ThreadPoolExecutor(max_workers=4) as executor:
        for symbol,model,error in executor.map(fetch,symbols):
            if model: models[symbol]=model
            else: failures[symbol]=error or 'canonical_coding_transcript_absent'
    variants,counts=[],Counter()
    for model in models.values():
        calls,stats=extract_gene(args.vcf,model);variants.extend(calls);counts.update(stats)
    ranked=rank_models(variants,prior)
    queried=vep_shortlist(cache,ranked,args.max_api_variants)
    ranked=rank_models(variants,prior)
    (work/'track1-candidates.json').write_text(json.dumps(ranked,indent=2))
    if not ranked: raise RuntimeError('No candidate model passed QC')
    winner=ranked[0];pair=winner['variants']
    for variant in pair:
        variant.setdefault('population_annotation','not_queried_frequency_unknown')
    sample=subprocess.check_output(['bcftools','query','-l',str(args.vcf)],text=True).strip().splitlines()
    if len(sample)!=1: raise ValueError('Exactly one sample required')
    columns=['proband_id','chrom_1','pos_1','ref_1','alt_1','chrom_2','pos_2','ref_2','alt_2','epcr','finding_type','notes']
    notes=f"{winner['gene']}; {winner['model']}; "+'; '.join(v['transcript']+':'+v['protein'] for v in pair)+'. Phase not established; parents unavailable. EPCR is subjective, not calibrated. Proband ID from VCF sample header.'
    row={'proband_id':sample[0],'epcr':args.epcr,'finding_type':'primary','notes':notes}
    for i,v in enumerate(pair,1): row.update({f'chrom_{i}':'chr'+v['chrom'],f'pos_{i}':v['pos'],f'ref_{i}':v['ref'],f'alt_{i}':v['alt']})
    with (output/f'track1_{args.team}.csv').open('w') as handle:
        writer=csv.DictWriter(handle,fieldnames=columns,lineterminator="\n");writer.writeheader();writer.writerow(row)
    nmd_results=[nmd(models[v['gene']],v) for v in pair if v['consequence']=='stop_gained']
    limits=['Phenotype-first shortlist, not exhaustive genome-wide consequence annotation.','Canonical coding transcripts only; structural variants and most splice variants not assessed.','Single-sample INFO/AF is not population frequency and is never used.','Rarity confirmed only where VEP returns frequencies; missing frequency is unknown, not zero.','Non-shortlisted variants use impact and phenotype heuristics without a rarity claim.','Recessive pairs unphased; dominant models merely de novo compatible without parental data.','All supplied terms included; one family-history term is withheld from the public copy.','EPCR is subjective, not a statistical output.']
    summary={'sample_id_source':'VCF sample header','assembly':'GRCh38','hpo_terms':terms,'hpo_information_weights':weights,'hpo_genes_ranked':len(prior),'gene_limit':args.gene_limit,'canonical_genes_scanned':len(models),'gene_model_failures':failures,'counts':dict(counts),'shortlist_api_variants':queried,'winning_model':winner,'winning_gene_hpo_rank':symbols.index(winner['gene'])+1,'hpo_matches':[t for t in terms if winner['gene'] in byterm[t]],'nmd':nmd_results,'limits':limits}
    (output/'track1_evidence.json').write_text(json.dumps(summary,indent=2))
    provenance=[]
    for meta in sorted((work/'track1-cache').glob('*.meta.json')):
        record=json.loads(meta.read_text())
        # Publish public-reference identities, not the private shortlist cache.
        if '/vep/' in record['url']:
            continue
        key=meta.name.split('.')[0]
        snapshot=work/'track1-cache'/(key+'.txt')
        if not snapshot.exists(): snapshot=work/'track1-cache'/(key+'.json')
        source=record['url'] if 'ensembl.org' in record['url'] else HPO_URL
        provenance.append({'source':source,'retrieved_utc':record['retrieved_utc'],'response_sha256':record['sha256'],'stored_snapshot_sha256':hashlib.sha256(snapshot.read_bytes()).hexdigest(),'bytes':record['bytes']})
    (output/'track1_reference_manifest.json').write_text(json.dumps(provenance,indent=2))
    (output/'track1_methods.md').write_text('# Phenotype-first discovery\n\nNo causal gene is supplied. Public HPO annotations rank genes by summed inverse prevalence of the eight terms; the top '+str(args.gene_limit)+' genes are scanned locally.\n\nQC: PASS, QUAL>=30, GQ>=30, DP>=10, biallelic calls; heterozygous balance 0.2–0.8. Canonical CDS references establish coding effect locally with reference matching and strand-aware translation. Score: HPO prior/3 + consequence (6 truncating; 2 missense) + inheritance compatibility (4 unphased recessive pair; 3 homozygous; 0 dominant). Up to 12 named shortlisted variants go to VEP; AF>1% excludes returned alleles.\n\n'+'\n'.join('- '+s for s in limits)+'\n\nCounts and NMD arithmetic: track1_evidence.json. NMD 50–55 nucleotide rule: PMID:35052820 (https://pubmed.ncbi.nlm.nih.gov/35052820/).\n')
    print(json.dumps({'winning_gene':winner['gene'],'model':winner['model'],'coding_candidates':len(variants),'genes_scanned':len(models),'api_variants':queried,'nmd':nmd_results},indent=2))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--vcf',type=Path,required=True);p.add_argument('--work',default='work');p.add_argument('--output',default='results');p.add_argument('--team',default='carson');p.add_argument('--hpo',action='append');p.add_argument('--gene-limit',type=int,default=40);p.add_argument('--max-api-variants',type=int,default=12,choices=range(1,13));p.add_argument('--epcr',type=float,default=0.90)
    args=p.parse_args()
    if not 0<args.epcr<=1: p.error('epcr must be in (0,1]')
    run(args)

if __name__=='__main__': main()
