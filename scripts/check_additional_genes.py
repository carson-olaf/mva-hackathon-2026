"""Reproduce the bounded ALPL/CDKN1C/HSD11B2 checks from local caches.

Default operation is offline. --fetch-approved permits only the ALPL/CDKN1C
canonical public annotations and requires explicit human download approval.
Individual records remain under ignored work/; stdout contains counts only.
This is a canonical-CDS check, not a clinical exclusion of these disorders.
"""
import argparse
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import requests
from mva_hackathon_2026.track1 import gene_model, extract_gene

ROOT = Path(__file__).resolve().parent.parent
VCF = ROOT / 'data/WGS_EX2312012_HGWCNDSX7.vcf.gz'
CACHE = ROOT / 'work/new-gene-public-cache'


class BoundedCache:
    def __init__(self, fetch):
        self.fetch = fetch
        self.total_bytes = 0
        self.allowed_transcripts = set()

    def get(self, url, params=None, text=False):
        assert not text
        key = hashlib.sha256(json.dumps([url, params], sort_keys=True).encode()).hexdigest()
        for base in [CACHE, ROOT/'work/track1-cache', ROOT/'work/clinical-cache']:
            path = base / (key + '.json')
            if path.exists():
                value = json.loads(path.read_text())
                self.register(value)
                return value
        allowed_gene = url in [f'https://rest.ensembl.org/lookup/symbol/homo_sapiens/{g}' for g in ['ALPL','CDKN1C']]
        allowed_cds = url in {f'https://rest.ensembl.org/sequence/id/{t}' for t in self.allowed_transcripts} and params == {'type':'cds'}
        if not self.fetch or not (allowed_gene or allowed_cds):
            raise RuntimeError('Unapproved or uncached reference: ' + url)
        with requests.get(url, params=params, headers={'Accept':'application/json','Content-Type':'application/json'}, stream=True, timeout=30) as response:
            response.raise_for_status()
            remaining = 300_000 - self.total_bytes
            if int(response.headers.get('Content-Length',0)) > remaining:
                raise RuntimeError('Approved public annotation transfer cap exceeded')
            chunks=[]
            for chunk in response.iter_content(4096):
                self.total_bytes += len(chunk)
                if self.total_bytes > 300_000:
                    raise RuntimeError('Approved public annotation transfer cap exceeded')
                chunks.append(chunk)
            raw=b''.join(chunks)
            value=json.loads(raw)
            CACHE.mkdir(exist_ok=True)
            (CACHE/(key+'.json')).write_bytes(raw)
            (CACHE/(key+'.meta.json')).write_text(json.dumps({'url':response.url,'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest(),'retrieved_utc':datetime.now(timezone.utc).isoformat()},indent=2))
        self.register(value)
        return value

    def register(self, value):
        if value.get('display_name') in ['ALPL','CDKN1C']:
            self.allowed_transcripts.update(t['id'] for t in value.get('Transcript',[]) if t.get('is_canonical') and t.get('Translation'))


def check(symbol, cache):
    destination = ROOT / f'work/new-{symbol.lower()}-check.json'
    if subprocess.run(['git', 'check-ignore', '-q', str(destination)], cwd=ROOT).returncode:
        raise RuntimeError('Patient intermediates must be git-ignored before this check runs')
    model=gene_model(cache,symbol)
    variants,counts=extract_gene(VCF,model)
    coding=set(model['positions'])
    boundaries={p for exon in model['exons'] for p in [exon['start'],exon['end']]}
    fmt='%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER[\t%GT\t%DP\t%GQ]\n'
    out=subprocess.run(['bcftools','query','-r',f"{model['chrom']}:{model['start']}-{model['end']}",'-f',fmt,str(VCF)],check=True,capture_output=True,text=True).stdout
    splice=[]
    for line in out.splitlines():
        chrom,pos,ref,alt,qual,filt,gt,dp,gq=line.split('\t')
        pos=int(pos)
        if pos in coding or not any(abs(pos-b)<=20 for b in boundaries): continue
        splice.append({'chrom':chrom,'pos':pos,'ref':ref,'alt':alt,'filter':filt,'gt':gt,'dp':dp,'gq':gq,'qual':qual,'distance_to_exon_edge':min(abs(pos-b) for b in boundaries)})
    result={'gene':symbol,'transcript':model['transcript'],'counts':dict(counts),'protein_altering_calls':variants,'noncoding_boundary_calls':splice,'limits':['Variant-only VCF; absence is not a reference genotype','Canonical CDS only; no SV, CNV, methylation, phase or biochemical assessment','Boundary records are not splice-effect predictions']}
    destination.write_text(json.dumps(result,indent=2))
    print(json.dumps({'gene':symbol,'transcript':model['transcript'],'counts':dict(counts),'protein_altering_call_count':len(variants),'noncoding_boundary_call_count':len(splice)},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--fetch-approved',action='store_true');p.add_argument('genes',nargs='+',choices=['ALPL','CDKN1C','HSD11B2']);args=p.parse_args()
    cache=BoundedCache(args.fetch_approved)
    for gene in args.genes: check(gene,cache)
    print(json.dumps({'public_annotation_bytes_fetched':cache.total_bytes}))
