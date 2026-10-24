"""Public-data BUB1B transcriptional reversal screen.

No patient data are read or transmitted. The null resamples complete unrelated
knockdown vectors within cell line, preserving gene covariance. Chemical
consensus profiles pool doses and cell lines. Provider checksum contradictions or
absent advertised checksums stop validated inference before correlations or nulls.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import json
import platform
from pathlib import Path
import re
import time

import numpy as np
import pandas as pd
import requests

BASE = "https://maayanlab.cloud/sigcom-lincs"
SHRNA = "8f1ff550-ece8-591d-a213-2763f854c008"
CRISPR = "96c7b8c5-1eca-5764-88e4-e4ccaee6603f"
MEAN_CP = "42cd56da-0ad8-5dad-b27c-fe1d135401b2"
CELLS = ("A375", "A549", "MCF7", "PC3")
MAX_BYTES = 1_000_000_000
OFFLINE = False


def cached(path: Path, url: str, payload: dict | None = None) -> Path:
    """Cache public bytes with an exact request and SHA-256 sidecar."""
    if path.exists() and path.stat().st_size:
        return path
    if OFFLINE:
        raise FileNotFoundError(f"Offline cache miss: {path.name}")
    path.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(4):
        try:
            response = requests.get(url, timeout=(20, 90), stream=True) if payload is None else requests.post(url, json=payload, timeout=(20, 90), stream=True)
            response.raise_for_status()
            if int(response.headers.get("Content-Length", "0")) > MAX_BYTES:
                raise ValueError("Download exceeds explicit 1 GB authorization boundary")
            size = 0
            digest = hashlib.sha256()
            temporary = path.with_suffix(path.suffix + ".part")
            with temporary.open("wb") as out:
                for chunk in response.iter_content(1024 * 1024):
                    size += len(chunk)
                    if size > MAX_BYTES:
                        raise ValueError("Streaming response exceeds 1 GB authorization boundary")
                    out.write(chunk)
                    digest.update(chunk)
            temporary.replace(path)
            path.with_name(path.name + ".source.json").write_text(json.dumps({"url": url, "payload": payload, "retrieved_utc": datetime.now(timezone.utc).isoformat(), "http_status": response.status_code, "content_type": response.headers.get("Content-Type"), "etag": response.headers.get("ETag"), "sha256": digest.hexdigest(), "bytes": size}, indent=2))
            return path
        except (requests.RequestException, OSError):
            if attempt == 3:
                raise
            time.sleep(attempt + 1)
    raise RuntimeError("Unreachable")


def metadata(cache: Path, label: str, where: dict, limit: int = 10000, skip: int = 0) -> list[dict]:
    path = cached(cache / f"{label}.json", BASE + "/metadata-api/signatures/find", {"filter": {"where": where, "limit": limit, "skip": skip}})
    return json.loads(path.read_text())


def normalize_name(name: str) -> str:
    return re.sub(r"[^a-z0-9]", "", name.casefold())


def approved_compounds(cache: Path) -> tuple[dict[str, dict], list[dict]]:
    pages = []
    for offset in range(0, 5000, 1000):
        path = cached(cache / f"chembl_phase4_page{offset}.json", f"https://www.ebi.ac.uk/chembl/api/data/molecule.json?max_phase=4&limit=1000&offset={offset}&only=molecule_chembl_id,pref_name,max_phase,molecule_synonyms")
        page = json.loads(path.read_text())
        pages.extend(page["molecules"])
        if offset + 1000 >= page["page_meta"]["total_count"]:
            break
    withdrawal = {}
    for offset in range(0, 5000, 1000):
        path = cached(cache / f"chembl_withdrawal{offset}.json", f"https://www.ebi.ac.uk/chembl/api/data/molecule.json?max_phase=4&limit=1000&offset={offset}&only=molecule_chembl_id,pref_name,withdrawn_flag,withdrawn_year,withdrawn_country,withdrawn_reason")
        page = json.loads(path.read_text())
        withdrawal.update({m["molecule_chembl_id"]: m.get("withdrawn_flag") for m in page["molecules"]})
        if offset + 1000 >= page["page_meta"]["total_count"]:
            break
    aliases = defaultdict(dict)
    for molecule in pages:
        if float(molecule.get("max_phase") or 0) != 4 or withdrawal.get(molecule["molecule_chembl_id"]) is not False:
            continue
        names = [molecule.get("pref_name") or ""] + [s.get("molecule_synonym") or s.get("synonyms") or "" for s in molecule.get("molecule_synonyms", [])]
        for name in names:
            key = normalize_name(name)
            if key:
                aliases[key][molecule["molecule_chembl_id"]] = molecule
    consensus = []
    for skip in range(0, 100000, 10000):
        # Earlier exploratory reads used these exact cache filenames.
        label = "mean_cp_metadata" if skip == 0 else f"mean_cp_metadata_{skip}"
        group = metadata(cache, label, {"library": MEAN_CP}, skip=skip)
        consensus.extend(group)
        if len(group) < 10000:
            break
    approved = {}
    for signature in consensus:
        name = signature["meta"]["pert_name"]
        matches = aliases.get(normalize_name(name), {})
        # Ambiguous normalized matches are rejected; salts are not guessed.
        if len(matches) == 1:
            molecule = next(iter(matches.values()))
            approved[name] = {"signature": signature, "chembl_id": molecule["molecule_chembl_id"], "approved_name": molecule["pref_name"], "max_phase": molecule["max_phase"]}
    return approved, consensus


def signature_path(cache: Path, signature: dict) -> Path:
    url = signature["meta"]["persistent_id"]
    suffix = ".tsv.gz" if url.endswith(".gz") else ".tsv"
    return cache / (signature["id"] + suffix)


def get_signature(cache: Path, signature: dict) -> Path:
    unavailable = cache / (signature["id"] + ".unavailable.json")
    if unavailable.exists():
        raise FileNotFoundError(f"Public signature unavailable; recorded in {unavailable.name}")
    try:
        return cached(signature_path(cache, signature), signature["meta"]["persistent_id"])
    except requests.HTTPError as error:
        if error.response is not None and error.response.status_code in (403, 404):
            unavailable.write_text(json.dumps({"id": signature["id"], "url": signature["meta"]["persistent_id"], "http_status": error.response.status_code}, indent=2))
        raise


def fetch_many(cache: Path, signatures: list[dict], workers: int) -> list[dict]:
    errors = []
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(get_signature, cache, sig): sig for sig in signatures}
        for count, future in enumerate(as_completed(futures), 1):
            sig = futures[future]
            try:
                future.result()
            except Exception as error:
                unavailable = cache / (sig["id"] + ".unavailable.json")
                errors.append(json.loads(unavailable.read_text()) if unavailable.exists() else {"id": sig["id"], "error": str(error)})
            if count % 100 == 0:
                print(f"Public signatures fetched {count}/{len(signatures)}; errors {len(errors)}", flush=True)
    return errors


def provider_integrity(cache: Path, signature: dict) -> dict:
    """Distinguish contradictory provenance from absent provider checksums."""
    path = signature_path(cache, signature)
    meta = signature["meta"]
    row = {"id": signature["id"], "library": signature["library"], "pert_name": meta.get("pert_name"), "url": meta["persistent_id"]}
    if not path.exists():
        return {**row, "status": "unavailable"}
    raw = path.read_bytes()
    actual_sha = hashlib.sha256(raw).hexdigest()
    actual_md5 = hashlib.md5(raw).hexdigest()
    advertised_sha = meta.get("sha256") or None
    advertised_md5 = meta.get("md5") or None
    advertised_size = meta.get("size_in_bytes") or None
    checks = []
    if advertised_sha:
        checks.append(advertised_sha == actual_sha)
    if advertised_md5:
        checks.append(advertised_md5 == actual_md5)
    if advertised_size:
        checks.append(advertised_size == len(raw))
    status = "mismatch" if checks and not all(checks) else "verified" if advertised_sha else "not_advertised"
    return {**row, "status": status, "actual_sha256": actual_sha, "advertised_sha256": advertised_sha, "actual_md5": actual_md5, "advertised_md5": advertised_md5, "actual_bytes": len(raw), "advertised_bytes": advertised_size}


def source_gate(cache: Path, signatures: list[dict]) -> dict:
    rows = [provider_integrity(cache, sig) for sig in signatures]
    counts = {}
    for row in rows:
        library_counts = counts.setdefault(row["library"], {})
        library_counts[row["status"]] = library_counts.get(row["status"], 0) + 1
    passed = bool(rows) and all(row["status"] == "verified" for row in rows)
    return {"evidence_status": "source_verified" if passed else "uninformative_quality_failure", "source_gate_passed": passed, "ranking_bonus": 0, "counts_by_library": counts, "reason": "Provider checksum contradictions or absent advertised checksums block validated cross-modality/drug inference. Local hashes and repeat HTTP equality prove transport/cache reproducibility, not identity with the advertised source version.", "records": rows}


def invalidate_previous_outputs(output: Path, provenance: dict) -> None:
    """Preserve old evidence while making its failed current validity explicit."""
    history = output / "signature_history"
    archived = []
    for name in ("signature_diagnostics.json", "signature_ranked_approved_drugs.csv", "signature_null_calibration.png"):
        path = output / name
        if not path.exists():
            continue
        raw = path.read_bytes()
        digest = hashlib.sha256(raw).hexdigest()
        history.mkdir(parents=True, exist_ok=True)
        saved = history / f"{path.stem}.{digest[:16]}{path.suffix}"
        if not saved.exists():
            saved.write_bytes(raw)
        archived.append({"original": name, "archive": saved.relative_to(output).as_posix(), "sha256": digest})
        invalid = {"evidence_status": "uninformative_quality_failure", "source_gate_passed": False, "current_valid": False, "ranking_bonus": 0, "ranking_eligible": False, "drug_reversal_inference_validated": False, "cross_modality_inference_validated": False, "historical_numeric_outputs_only": True}
        if path.suffix == ".json":
            try:
                previous = json.loads(raw)
            except json.JSONDecodeError:
                previous = {}
            if not isinstance(previous, dict):
                previous = {}
            previous.update(invalid)
            previous["invalidation_reason"] = provenance["reason"]
            for row in previous.get("top10", []):
                if isinstance(row, dict):
                    row.update(invalid)
            path.write_text(json.dumps(previous, indent=2))
        elif path.suffix == ".csv":
            try:
                previous_table = pd.read_csv(path)
            except (pd.errors.EmptyDataError, pd.errors.ParserError):
                # The original remains in history, without a misleading current
                # filename. The run status links to the preserved artifact.
                path.unlink()
                continue
            for field, value in invalid.items():
                previous_table[field] = value
            previous_table.to_csv(path, index=False)
        else:
            # An image cannot carry machine-readable validity fields. Remove
            # its current filename after preserving its original bytes.
            path.unlink()
    status = {"evidence_status": "uninformative_quality_failure", "source_gate_passed": False, "current_valid": False, "ranking_bonus": 0, "exit_code": 2, "reason": provenance["reason"], "historical_artifacts": archived}
    (output / "signature_run_status.json").write_text(json.dumps(status, indent=2))


def read_vector(path: Path, genes: list[str]) -> np.ndarray:
    data = pd.read_csv(path, sep="\t", index_col=0)
    if data.shape[1] != 1:
        raise ValueError(f"Expected single public expression vector: {path.name}")
    series = pd.to_numeric(data.iloc[:, 0], errors="raise")
    if series.index.has_duplicates:
        series = series.groupby(level=0).mean()
    vector = series.reindex(genes).to_numpy(dtype=float)
    if not np.isfinite(vector).all():
        raise ValueError(f"Missing landmark expression genes in {path.name}")
    return vector


def centered_rank(vector: np.ndarray) -> np.ndarray:
    ranks = pd.Series(vector).rank(method="average").to_numpy(dtype=float)
    ranks = ranks - ranks.mean()
    norm = np.linalg.norm(ranks)
    if not np.isfinite(norm) or norm == 0:
        raise ValueError("Cannot correlate a non-finite or constant expression vector")
    return ranks / norm


def cell_consensus(cache: Path, signatures: list[dict], genes: list[str]) -> dict[str, np.ndarray]:
    vectors = defaultdict(list)
    for sig in signatures:
        cell = sig["meta"].get("cell_line")
        if cell not in CELLS or not signature_path(cache, sig).exists():
            continue
        vectors[cell].append(centered_rank(read_vector(signature_path(cache, sig), genes)))
    return {cell: np.mean(rows, axis=0) for cell, rows in vectors.items()}


def bh(p: np.ndarray) -> np.ndarray:
    order = np.argsort(p)
    adjusted = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    result = np.empty_like(p)
    result[order] = np.minimum(adjusted, 1)
    return result


def verified_shrna_reproducibility(cache: Path, output: Path, resamples: int = 999, seed: int = 20261024) -> list[dict]:
    """Assess within-cell reproducibility using only provider-verified vectors.

    This diagnostic does not consume drug or CRISPR signatures. Whole profiles
    are resampled within cell/time to preserve covariance between landmark genes.
    """
    gene_file = cached(cache / "geneinfo_beta.txt", "https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/geneinfo_beta.txt")
    gd = pd.read_csv(gene_file, sep="\t")
    genes = sorted(set(gd.loc[gd.feature_space == "landmark", "gene_symbol"]) - {"BUB1B"})
    bub = metadata(cache, "bub1b_metadata", {"meta.pert_name": "BUB1B"}, limit=200)
    control_rng = np.random.default_rng(seed)
    targets = sorted(control_rng.choice(genes, size=100, replace=False).tolist())
    controls = metadata(cache, f"controls_100_{seed}", {"library": SHRNA, "meta.pert_name": {"inq": targets}, "meta.cell_line": {"inq": list(CELLS)}}, limit=10000)
    rng = np.random.default_rng(seed)
    rows = []
    input_ids = {}
    for cell in CELLS:
        targets_in_cell = [x for x in bub if x["library"] == SHRNA and x["meta"].get("pert_time") == "96 h" and x["meta"].get("cell_line") == cell]
        control_profiles = [x for x in controls if x["meta"].get("pert_time") == "96 h" and x["meta"].get("cell_line") == cell and signature_path(cache, x).exists()]
        def vector(sig: dict) -> np.ndarray:
            audit = provider_integrity(cache, sig)
            if audit["status"] != "verified":
                raise ValueError(f"Provider integrity failed in shRNA diagnostic: {sig['id']}")
            return centered_rank(read_vector(signature_path(cache, sig), genes))
        target_matrix = np.stack([vector(sig) for sig in targets_in_cell])
        pool = np.stack([vector(sig) for sig in control_profiles])
        n = len(target_matrix)
        triangle = np.triu_indices(n, 1)
        observed = float(np.median((target_matrix @ target_matrix.T)[triangle]))
        null = []
        for _ in range(resamples):
            sampled = pool[rng.choice(len(pool), n, replace=False)]
            null.append(float(np.median((sampled @ sampled.T)[triangle])))
        rows.append({"cell_line": cell, "BUB1B_profiles": n, "BUB1B_pairs": len(triangle[0]), "landmarks_excluding_BUB1B": len(genes), "observed_median_pairwise_spearman": observed, "null_resamples": resamples, "same_cell_unrelated_shRNA_profile_pool": len(pool), "null_median": float(np.median(null)), "null_95th_percentile": float(np.quantile(null, .95)), "competitive_reproducibility_p": float((1 + np.sum(np.array(null) >= observed)) / (resamples + 1)), "seed": seed, "provider_sha256_and_size_all_passed": True})
        input_ids[cell] = {"BUB1B_signature_ids": [x["id"] for x in targets_in_cell], "null_pool_signature_ids": [x["id"] for x in control_profiles]}
    table = pd.DataFrame(rows)
    table["BH_q_across_four_cells"] = bh(table.competitive_reproducibility_p.to_numpy())
    output.mkdir(parents=True, exist_ok=True)
    table.to_csv(output / "signature_shRNA_verified_reproducibility.csv", index=False)
    (output / "signature_shRNA_qc_input_ids.json").write_text(json.dumps(input_ids, indent=2))
    print(table.to_string(index=False), flush=True)
    return table.to_dict(orient="records")


def run(cache: Path, output: Path, permutations: int = 4999, workers: int = 6, control_genes: int = 100, seed: int = 20261024) -> bool:
    if permutations < 99 or not 1 <= workers <= 12 or control_genes < 20:
        raise ValueError("Require >=99 resamples, 1-12 workers, and >=20 unrelated control genes")
    cache.mkdir(parents=True, exist_ok=True)
    output.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    gene_file = cached(cache / "geneinfo_beta.txt", "https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/geneinfo_beta.txt")
    gene_data = pd.read_csv(gene_file, sep="\t")
    genes = sorted(set(gene_data.loc[gene_data.feature_space == "landmark", "gene_symbol"]) - {"BUB1B"})
    bub = metadata(cache, "bub1b_metadata", {"meta.pert_name": "BUB1B"}, limit=200)
    shrna = [s for s in bub if s["library"] == SHRNA and s["meta"].get("cell_line") in CELLS and s["meta"].get("pert_time") == "96 h"]
    crispr = [s for s in bub if s["library"] == CRISPR and s["meta"].get("cell_line") in CELLS and s["meta"].get("pert_time") == "96 h"]
    candidate_controls = sorted(rng.choice(genes, size=min(control_genes, len(genes)), replace=False).tolist())
    controls = metadata(cache, f"controls_{control_genes}_{seed}", {"library": SHRNA, "meta.pert_name": {"inq": candidate_controls}, "meta.cell_line": {"inq": list(CELLS)}}, limit=10000)
    controls = [s for s in controls if s["meta"].get("pert_time") == "96 h"]
    approved, all_compounds = approved_compounds(cache)
    print(f"BUB1B shRNA {len(shrna)}, CRISPR {len(crispr)}; control profiles {len(controls)}; approved consensus drugs {len(approved)}/{len(all_compounds)}", flush=True)
    # Limit each unrelated target/cell to three hairpin profiles, chosen by ID,
    # rather than weight targets with more published replicates more heavily.
    selected = []
    grouped = defaultdict(list)
    for sig in controls:
        grouped[(sig["meta"]["pert_name"], sig["meta"]["cell_line"])].append(sig)
    for group in grouped.values():
        selected.extend(sorted(group, key=lambda s: s["id"])[:3])
    errors = fetch_many(cache, shrna + crispr + selected + [x["signature"] for x in approved.values()], workers)
    (cache / "fetch_errors.json").write_text(json.dumps(errors, indent=2))
    provenance = source_gate(cache, shrna + crispr + selected + [x["signature"] for x in approved.values()])
    (output / "signature_source_qc.json").write_text(json.dumps(provenance, indent=2))
    if not provenance["source_gate_passed"]:
        # Invalidate any earlier success before returning, including artifacts
        # produced by a prior run with a different source version.
        invalidate_previous_outputs(output, provenance)
        print(json.dumps({k: v for k, v in provenance.items() if k != "records"}, indent=2), flush=True)
        return False
    sh_cells = cell_consensus(cache, shrna, genes)
    cr_cells = cell_consensus(cache, crispr, genes)
    common_cells = [c for c in CELLS if c in sh_cells and c in cr_cells]
    if len(common_cells) < 3:
        raise ValueError("Fewer than three matched cells with BUB1B shRNA and CRISPR")
    control_cells = defaultdict(list)
    coherent_controls = defaultdict(dict)
    for (target, cell), group in grouped.items():
        available = [s for s in sorted(group, key=lambda s: s["id"])[:3] if signature_path(cache, s).exists()]
        if available and cell in common_cells:
            value = np.mean([centered_rank(read_vector(signature_path(cache, s), genes)) for s in available], axis=0)
            control_cells[cell].append(value)
            coherent_controls[target][cell] = value
    sh = centered_rank(np.mean([sh_cells[c] for c in common_cells], axis=0))
    cr = centered_rank(np.mean([cr_cells[c] for c in common_cells], axis=0))
    drugs, vectors = [], []
    for name, info in sorted(approved.items()):
        path = signature_path(cache, info["signature"])
        if not path.exists():
            continue
        try:
            vector = centered_rank(read_vector(path, genes))
        except ValueError as error:
            errors.append({"drug": name, "error": str(error)})
            continue
        drugs.append({"drug": name, **{k: info[k] for k in ("chembl_id", "approved_name", "max_phase")}, "signature_id": info["signature"]["id"], "signature_url": info["signature"]["meta"]["persistent_id"], "withdrawn_flag": False})
        vectors.append(vector)
    chemical = np.stack(vectors).astype(np.float64)
    # Positive reversal is anti-correlation. Require both genetic modalities to
    # reverse: the minimum prevents one modality driving the discovery statistic.
    sh_score, cr_score = -(chemical @ sh), -(chemical @ cr)
    discovery_score = np.minimum(sh_score, cr_score)
    count = len(drugs)
    exceed_sh = np.zeros(count)
    exceed_min = np.zeros(count)
    exceed_family = np.zeros(count)
    maxima = []
    # Cell-stratified target-label randomization samples whole knockdown vectors.
    # It retains correlations between landmark genes within a cell profile. It
    # is a competitive gene-specificity null, not proof of treatment response.
    for start in range(0, permutations, 100):
        n = min(100, permutations - start)
        null_vectors = np.zeros((n, len(genes)))
        for cell in common_cells:
            pool = np.stack(control_cells[cell])
            null_vectors += pool[rng.integers(len(pool), size=n)] / len(common_cells)
        null_vectors = np.stack([centered_rank(x) for x in null_vectors])
        null_scores = -(chemical @ null_vectors.T)
        exceed_sh += (null_scores >= sh_score[:, None]).sum(axis=1)
        # Conservative comparison: random shRNA signature must exceed the
        # smaller observed shRNA/CRISPR reversal score.
        exceed_min += (null_scores >= discovery_score[:, None]).sum(axis=1)
        null_max = null_scores.max(axis=0)
        exceed_family += (null_max[None, :] >= discovery_score[:, None]).sum(axis=1)
        maxima.extend(null_max.tolist())
    p_sh = (exceed_sh + 1) / (permutations + 1)
    p_convergence = (exceed_min + 1) / (permutations + 1)
    p_family = (exceed_family + 1) / (permutations + 1)
    result = pd.DataFrame(drugs)
    result["shRNA_reversal_rho"] = sh_score
    result["CRISPR_reversal_rho"] = cr_score
    result["minimum_modality_reversal_rho"] = discovery_score
    result["competitive_shRNA_p"] = p_sh
    result["competitive_conservative_p"] = p_convergence
    result["bh_q"] = bh(p_convergence)
    result["maxT_familywise_p"] = p_family
    result = result.sort_values("minimum_modality_reversal_rho", ascending=False)
    result.insert(0, "rank", np.arange(1, len(result) + 1))
    result.to_csv(output / "signature_ranked_approved_drugs.csv", index=False)
    # Coherent unrelated targets are a sensitivity check on the cell-wise null.
    coherent = [centered_rank(np.mean(list(cells.values()), axis=0)) for cells in coherent_controls.values() if set(common_cells).issubset(cells)]
    coherent_scores = -(chemical @ np.stack(coherent).T) if coherent else np.empty((count, 0))
    negative_target_profiles = {"shRNA": [], "CRISPR": []}
    target_qc = {}
    for name, profiles in (("shRNA", shrna), ("CRISPR", crispr)):
        values = []
        for sig in profiles:
            path = signature_path(cache, sig)
            if path.exists():
                data = pd.read_csv(path, sep="\t", index_col=0).iloc[:, 0]
                if "BUB1B" in data.index:
                    value = float(data.loc["BUB1B"])
                    values.append(value)
                    if value < 0:
                        negative_target_profiles[name].append(sig)
        target_qc[name] = {"profiles": len(values), "negative_target_coefficient": sum(v < 0 for v in values), "median_target_coefficient": float(np.median(values)), "target_is_measured_landmark": True, "caveat": "A signed characteristic-direction coefficient is not a fold-change measurement or direct knockdown-efficiency assay."}
    negative_sh_cells = cell_consensus(cache, negative_target_profiles["shRNA"], genes)
    negative_cr_cells = cell_consensus(cache, negative_target_profiles["CRISPR"], genes)
    sensitivity_cells = [c for c in common_cells if c in negative_sh_cells and c in negative_cr_cells]
    negative_target_rho = float(centered_rank(np.mean([negative_sh_cells[c] for c in sensitivity_cells], axis=0)) @ centered_rank(np.mean([negative_cr_cells[c] for c in sensitivity_cells], axis=0)))
    diagnostics = {
        "seed": seed, "permutations": permutations, "public_only": True,
        "software_versions": {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__, "requests": requests.__version__},
        "analysis_code_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "evidence_status": "source_verified",
        "source_gate_passed": True,
        "target_expression_QC": target_qc,
        "sensitivity_negative_target_shRNA_CRISPR_rho": negative_target_rho,
        "cell_lines": common_cells, "perturbation_time": "96 h", "landmark_genes": len(genes),
        "BUB1B_shRNA_profiles": len(shrna), "BUB1B_CRISPR_profiles": len(crispr),
        "BUB1B_shRNA_CRISPR_rho": float(sh @ cr),
        "same_cell_shRNA_CRISPR_rho": {c: float(centered_rank(sh_cells[c]) @ centered_rank(cr_cells[c])) for c in common_cells},
        "control_gene_targets_requested": control_genes,
        "control_target_cell_counts": {c: len(control_cells[c]) for c in common_cells},
        "coherent_control_targets_all_cells": len(coherent),
        "chemical_consensus_profiles_metadata": len(all_compounds),
        "approved_exact_unique_ChEMBL_matches": len(approved),
        "approved_drugs_analyzed": len(result), "fetch_errors": len(errors),
        "positive_both_modalities": int((discovery_score > 0).sum()),
        "BH_q_below_0_05": int((result.bh_q < .05).sum()),
        "maxT_familywise_p_below_0_05": int((result.maxT_familywise_p < .05).sum()),
        "null_max_95pct": float(np.quantile(maxima, .95)),
        "top10": result.head(10).to_dict(orient="records"),
        "shRNA_cell_pair_rho": {f"{a}:{b}": float(centered_rank(sh_cells[a]) @ centered_rank(sh_cells[b])) for i,a in enumerate(common_cells) for b in common_cells[i+1:]},
        "CRISPR_cell_pair_rho": {f"{a}:{b}": float(centered_rank(cr_cells[a]) @ centered_rank(cr_cells[b])) for i,a in enumerate(common_cells) for b in common_cells[i+1:]},
    }
    # Preserve specificity against coherent gene controls for top-ranked drugs.
    for row in diagnostics["top10"]:
        index = next(i for i, drug in enumerate(drugs) if drug["drug"] == row["drug"])
        row["coherent_control_competitive_p"] = float((1 + np.sum(coherent_scores[index] >= discovery_score[index])) / (len(coherent) + 1))
    (output / "signature_diagnostics.json").write_text(json.dumps(diagnostics, indent=2))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    axes[0].scatter(sh_score, cr_score, s=9, alpha=.45, color="#164d7b")
    axes[0].axhline(0, color="black", lw=.7)
    axes[0].axvline(0, color="black", lw=.7)
    axes[0].set(xlabel="shRNA reversal (negative Spearman rho)", ylabel="CRISPR reversal (negative Spearman rho)", title=f"{len(result):,} approved-drug consensus profiles")
    axes[1].hist(maxima, bins=35, color="#164d7b", alpha=.8)
    axes[1].axvline(float(result.iloc[0].minimum_modality_reversal_rho), color="#b43d32", lw=2, label="Best observed two-modality score")
    axes[1].set(xlabel="Maximum reversal across drugs", ylabel="Null resamples", title="Cell-stratified unrelated-knockdown null")
    axes[1].legend(fontsize=8)
    fig.suptitle("Exploratory reversal does not demonstrate correction of chromosome instability", fontsize=11)
    fig.tight_layout()
    fig.savefig(output / "signature_null_calibration.png", dpi=180)
    plt.close(fig)
    manifest = []
    for file in sorted(cache.iterdir()):
        if file.is_file() and not file.name.endswith((".source.json", ".part", ".log")):
            manifest.append({"file": file.name, "sha256": hashlib.sha256(file.read_bytes()).hexdigest(), "bytes": file.stat().st_size})
    (output / "signature_public_cache_manifest.json").write_text(json.dumps(manifest, indent=2))
    print(json.dumps({k:v for k,v in diagnostics.items() if k not in ("top10", "shRNA_cell_pair_rho", "CRISPR_cell_pair_rho")}, indent=2), flush=True)
    print(result.head(10).to_string(index=False), flush=True)
    (output / "signature_run_status.json").write_text(json.dumps({"evidence_status": "source_verified", "source_gate_passed": True, "current_valid": True, "exit_code": 0}, indent=2))
    return True


def main() -> None:
    global OFFLINE
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=Path("work/signature"))
    parser.add_argument("--output", type=Path, default=Path("results"))
    parser.add_argument("--permutations", type=int, default=4999)
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--control-genes", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20261024)
    parser.add_argument("--offline", action="store_true", help="Require cached public responses and never use network")
    parser.add_argument("--shrna-qc", action="store_true", help="Run only provider-verified shRNA reproducibility with 999 whole-profile null resamples")
    args = parser.parse_args()
    OFFLINE = args.offline
    if args.shrna_qc:
        verified_shrna_reproducibility(args.cache, args.output, seed=args.seed)
    else:
        passed = run(args.cache, args.output, args.permutations, args.workers, args.control_genes, args.seed)
        if not passed:
            raise SystemExit(2)


if __name__ == "__main__":
    main()
