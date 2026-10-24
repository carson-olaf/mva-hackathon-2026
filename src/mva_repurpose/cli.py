import argparse
import hashlib
import json
from pathlib import Path

from .cache import PublicCache
from .engine import load_pgx, run, write_outputs


def main(argv=None):
    parser = argparse.ArgumentParser(description="Direction-aware rare-disease drug hypotheses from public APIs only")
    parser.add_argument("--gene", help="HGNC symbol; infer from named HGVS variants if omitted")
    parser.add_argument("--variant", action="append", default=[], help="Named HGVS variant (at most ten); repeat as needed")
    parser.add_argument("--hpo", action="append", default=[], help="HP:nnnnnnn; repeat as needed")
    parser.add_argument("--mechanism", choices=["loss_of_function", "gain_of_function", "unknown"], default="unknown")
    parser.add_argument("--disease", help="Exact disease name for matching ChEMBL phase-4 indication labels")
    parser.add_argument("--pgx", type=Path, help="Small gene/phenotype/candidate_drugs_affected table, CSV/TSV or JSON")
    parser.add_argument("--evidence", type=Path, help="Explicit curated literature bridges; see examples/evidence-schema.json")
    parser.add_argument("--approval-audit", type=Path, help="Independent official-product approval verification JSON for curated hypotheses")
    parser.add_argument("--cache", type=Path, default=Path("work/repurpose-cache"))
    parser.add_argument("--offline", action="store_true", help="Forbid all network calls; fail on any cache miss")
    parser.add_argument("--cache-fallback", action="append", default=[], type=Path, help="Read-only public cache to reuse without duplicating references")
    parser.add_argument("--output", type=Path, default=Path("results/repurpose"), help="Output prefix for JSON, CSV and Markdown")
    parser.add_argument("--max-pathways", type=int, default=20)
    parser.add_argument("--max-proteins", type=int, default=200)
    args = parser.parse_args(argv)
    if not 1 <= args.max_pathways <= 100 or not 1 <= args.max_proteins <= 500:
        parser.error("Pathway cap must be 1–100; protein cap must be 1–500")
    try:
        evidence = json.loads(args.evidence.read_text()) if args.evidence else None
        approval = json.loads(args.approval_audit.read_text()) if args.approval_audit else None
        result = run(PublicCache(args.cache, offline=args.offline, fallbacks=tuple(args.cache_fallback)), gene=args.gene,
                     variants=args.variant, hpo=args.hpo, mechanism=args.mechanism,
                     disease=args.disease, pgx=load_pgx(args.pgx), evidence=evidence, approval_audit=approval,
                     max_pathways=args.max_pathways, max_proteins=args.max_proteins)
        result["optional_input_provenance"] = {
            name: {"file": path.name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
            for name, path in (("curated_evidence", args.evidence), ("approval_audit", args.approval_audit), ("pgx", args.pgx))
            if path is not None
        }
        result["software_provenance"] = {
            source.name: hashlib.sha256(source.read_bytes()).hexdigest()
            for source in sorted(Path(__file__).parent.glob("*.py"))
        }
        write_outputs(result, args.output)
    except (ValueError, FileNotFoundError) as exc:
        parser.error(str(exc))
    print(json.dumps({"gene": result["gene"], "phase4_proxy_graph_candidates": len(result["ranking"]),
                      "supported_hypotheses": len(result["supported_hypotheses"]),
                      "known_disease_label_matches": result["known_disease_label_matches"],
                      "pgx_matched_rows": result["pgx_integration"]["matched_input_rows"],
                      "curated_laboratory_hypotheses": len(result["curated_hypotheses"]["ranked"]),
                      "abstention": result["abstention"], "output": str(args.output)}))
    return 0
