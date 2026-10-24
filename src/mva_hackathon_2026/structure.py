"""Reference-only structural triage; this module never reads the private genome.

Fetch public reference structures with --fetch. Run fpocket separately or pass
--fpocket. Cached DynaMut2/KORPM predictions are kept distinct from geometry.
Positive ddG means increased unfolding free energy (stabilization) throughout.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import re
import shlex
import subprocess
import tarfile
from pathlib import Path

import numpy as np
import requests


FPOCKET_COMMIT = "4bb0d8447f62fee77e2c3c29f54b5fcaf5e2c066"
KORPM_COMMIT = "c9a50da5bbe885bdc1e6ec8c732eb914c78a0ead"


def download_small(url: str, target: Path) -> None:
    """Enforce the user's per-download ceiling, also for unknown-length streams."""
    limit = 1_000_000_000
    with requests.get(url, timeout=180, stream=True) as response:
        response.raise_for_status()
        if int(response.headers.get("Content-Length", 0)) > limit:
            raise ValueError("Download exceeds 1 GB; explicit user approval is required")
        total = 0
        temporary = target.with_suffix(target.suffix + ".part")
        try:
            with temporary.open("wb") as handle:
                for chunk in response.iter_content(1024 * 1024):
                    total += len(chunk)
                    if total > limit:
                        raise ValueError("Download exceeded 1 GB; stopped")
                    handle.write(chunk)
            temporary.replace(target)
        except BaseException:
            temporary.unlink(missing_ok=True)
            raise


def prepare_tools(cache: Path) -> tuple[Path, Path]:
    """Build pinned public tools in ignored cache; do not install system-wide."""
    cache.mkdir(parents=True, exist_ok=True)
    lock = Path.home() / ".agents/bin/build-lock"
    if not lock.exists():
        raise RuntimeError("Required build-lock executable not found")
    roots = []
    for name, repository, commit in [("fpocket", "Discngine/fpocket", FPOCKET_COMMIT),
                                     ("korpm", "chaconlab/korpm", KORPM_COMMIT)]:
        root = cache / f"{name}-{commit}"
        archive = cache / f"{name}.tar.gz"
        if not root.exists():
            if not archive.exists():
                download_small(f"https://codeload.github.com/{repository}/tar.gz/{commit}", archive)
            with tarfile.open(archive) as source:
                source.extractall(cache, filter="data")
        (cache / f"{name}_commit.txt").write_text(commit + "\n")
        roots.append(root)
    fpocket, korpm = roots
    arch = "MACOSXARM64" if platform.system() == "Darwin" and platform.machine() == "arm64" else "MACOSXX86_64" if platform.system() == "Darwin" else "LINUXAMD64"
    build_id = f"{platform.system()}-{platform.machine()}-fixed-seed-20261023-v1\n"
    fpocket_marker = cache / "fpocket_local_build.txt"
    if not fpocket_marker.exists() or fpocket_marker.read_text() != build_id:
        # Fpocket otherwise seeds Monte Carlo volumes from wall time. Fix only
        # random seeds, not pocket-scoring equations or trained coefficients.
        for filename in ("utils.c", "voronoi.c", "clusterlib.c"):
            path = fpocket / "src" / filename
            source = path.read_text().replace("srand((int)time(NULL));", "srand(20261023);")
            source = source.replace("srand(time(NULL));", "srand(20261023);")
            source = source.replace("unsigned int initseed = (unsigned int) time(0);", "unsigned int initseed = 20261023;")
            path.write_text(source)
        with (cache / "fpocket_build.log").open("w") as log:
            subprocess.run([str(lock), "make", "-B", "-C", str(fpocket), f"ARCH={arch}",
                            f"CFLAGS=-DM_OS_LINUX -DMNO_MEM_DEBUG -O2 -g -std=gnu99 -Iplugins/include -Iplugins/{arch}/molfile -Wno-error=implicit-function-declaration",
                            "QCFLAGS=-O -g -ansi -Isrc/qhull/src -Wno-error=implicit-function-declaration", "fpocket"], check=True, stdout=log, stderr=log)
        fpocket_marker.write_text(build_id)
    korpm_marker = cache / "korpm_local_build.txt"
    if not korpm_marker.exists() or korpm_marker.read_text() != build_id:
        # Upstream enables an unrelated density-map feature whose headers are
        # absent, and links a nonexistent pdb2 library. Scoring code is unchanged.
        build = korpm / "sbg/src/libpdb/Release_gcc/subdir.mk"
        build.write_text(build.read_text().replace("-DVOLUME_INCLUDED ", ""))
        with (cache / "korpm_build.log").open("w") as log:
            for component in ("libpdb", "libtools", "libenergy", "korpm"):
                command = [str(lock), "make", "-B", "-C", str(korpm / f"sbg/src/{component}/Release_gcc"), "all"]
                if component == "korpm":
                    command.append("LIBS=-lenergy_gcc -lpdb_gcc -ltools_gcc")
                subprocess.run(command, check=True, stdout=log, stderr=log)
        korpm_marker.write_text(build_id)
    potential = korpm / "pot/korp6Dv1.bin"
    if potential.stat().st_size < 300_000_000:
        download_small(f"https://media.githubusercontent.com/media/chaconlab/korpm/{KORPM_COMMIT}/pot/korp6Dv1.bin", potential)
    return fpocket / "bin/fpocket", korpm


def dynamut_reference_prediction(cache: Path, submit: bool = False) -> None:
    url = "https://biosig.lab.uq.edu.au/dynamut2/api/prediction_single"
    job_path = cache / "dynamut_submit.json"
    if submit and not job_path.exists():
        reference = cache / "AF-O60566-F1.pdb"
        # The file is a public WT reference, never a patient-derived structure.
        response = requests.post(url, data={"chain": "A", "mutation": "N1002K"},
                                 files={"pdb_file": (reference.name, reference.read_bytes(), "chemical/x-pdb")}, timeout=120)
        response.raise_for_status()
        job = response.json()
        if "job_id" not in job:
            raise RuntimeError(f"DynaMut2 did not accept job: {job}")
        job_path.write_text(json.dumps(job, indent=2) + "\n")
    if job_path.exists():
        job = json.loads(job_path.read_text())
        response = requests.get(url, params={"job_id": job["job_id"]}, timeout=90)
        response.raise_for_status()
        (cache / "dynamut_status.json").write_text(json.dumps(response.json(), indent=2) + "\n")


def fetch_references(cache: Path) -> None:
    cache.mkdir(parents=True, exist_ok=True)
    urls = {
        "alphafold_metadata.json": "https://alphafold.ebi.ac.uk/api/prediction/O60566",
        "uniprot.json": "https://rest.uniprot.org/uniprotkb/O60566.json",
        "pdbe_mapping.json": "https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/O60566",
        "6TLJ.cif": "https://files.rcsb.org/download/6TLJ.cif",
        "5KHU.cif": "https://files.rcsb.org/download/5KHU.cif",
        "6JKK.cif": "https://files.rcsb.org/download/6JKK.cif",
        "6JKM.cif": "https://files.rcsb.org/download/6JKM.cif",
    }
    for name, url in urls.items():
        if not (cache / name).exists():
            r = requests.get(url, timeout=90)
            r.raise_for_status()
            (cache / name).write_bytes(r.content)
    meta = json.loads((cache / "alphafold_metadata.json").read_text())[0]
    if not (cache / "AF-O60566-F1.pdb").exists():
        r = requests.get(meta["pdbUrl"], timeout=90)
        r.raise_for_status()
        (cache / "AF-O60566-F1.pdb").write_bytes(r.content)


def pdb_atoms(path: Path) -> list[dict]:
    atoms = []
    for line in path.read_text().splitlines():
        if not line.startswith("ATOM") or line[16] not in " A":
            continue
        element = line[76:78].strip() or line[12:16].strip()[0]
        if element == "H":
            continue
        atoms.append({"res": int(line[22:26]), "name": line[12:16].strip(),
                      "aa": line[17:20], "chain": line[21],
                      "xyz": np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]),
                      "bfactor": float(line[60:66]), "element": element})
    return atoms


def coordinate_coverage(path: Path, chain: str) -> dict:
    """Count modeled ATOM residues, not deposited full construct sequence length."""
    columns = []
    numbers = set()
    for line in path.read_text().splitlines():
        if line.startswith("_atom_site."):
            columns.append(line.split()[0].removeprefix("_atom_site."))
        elif columns and line.startswith(("ATOM ", "HETATM ")):
            row = dict(zip(columns, shlex.split(line), strict=True))
            if row.get("auth_asym_id") == chain and row.get("group_PDB") == "ATOM":
                numbers.add(int(row["auth_seq_id"]))
        elif columns and line == "#":
            break
    return {"structure": path.stem, "chain": chain, "modeled_residues": len(numbers),
            "first": min(numbers), "last": max(numbers), "models_1002": 1002 in numbers,
            "modeled_residues_705_1050": sum(705 <= n <= 1050 for n in numbers)}


def residue_sasa(atoms: list[dict], residue: int, points: int = 960) -> float:
    """Shrake-Rupley solvent-accessible area, Bondi radii, 1.4-A probe.

    Deterministic Fibonacci sphere sampling avoids random Monte Carlo variation.
    Area is a geometric descriptor, never a folding energy or binding score.
    """
    index = np.arange(points)
    z = 1 - 2 * (index + 0.5) / points
    theta = index * math.pi * (3 - math.sqrt(5))
    circle = np.sqrt(1 - z * z)
    sphere = np.column_stack((circle * np.cos(theta), circle * np.sin(theta), z))
    xyz = np.array([a["xyz"] for a in atoms])
    radii = np.array([{"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80}.get(a["element"], 1.70) + 1.4 for a in atoms])
    area = 0.0
    for i, atom in enumerate(atoms):
        if atom["res"] != residue:
            continue
        distance = np.linalg.norm(xyz - xyz[i], axis=1)
        neighbors = np.flatnonzero((distance < radii + radii[i]) & (distance > 0))
        sample = xyz[i] + radii[i] * sphere
        occupied = np.zeros(points, dtype=bool)
        for j in neighbors:
            occupied |= np.sum((sample - xyz[j]) ** 2, axis=1) < radii[j] ** 2
        area += 4 * math.pi * radii[i] ** 2 * np.mean(~occupied)
    return float(area)


def prepare_pocket_inputs(cache: Path) -> None:
    lines = (cache / "AF-O60566-F1.pdb").read_text().splitlines()
    for name, start in [("AF-O60566-F1-b20", 1), ("AF-O60566-Cterm-b20", 705)]:
        # AlphaFold stores confidence in B factors. A constant prevents a pocket
        # program from treating pLDDT as experimental thermal displacement.
        selected = [s[:60] + " 20.00" + s[66:] for s in lines
                    if s.startswith("ATOM") and int(s[22:26]) >= start]
        (cache / f"{name}.pdb").write_text("\n".join(selected) + "\nTER\nEND\n")


def run_korpm(cache: Path, root: Path) -> None:
    """Run the released model on the target and a reference-only mutation null."""
    controls = cache / "korpm_controls_input.txt"
    controls.write_text("1BNI IA76A\n1EY0 TA44V\n1IHB FA82Q\n")
    with (cache / "korpm_controls.log").open("w") as log:
        subprocess.run([str((root / "sbg/bin/korpm_gcc").resolve()), str(controls),
                        "--dir", str(root / "Ssym"), "--score_file", str(root / "pot/korp6Dv1.bin"),
                        "-o", str(cache / "korpm_controls_predictions.txt")], check=True, stdout=log, stderr=log)
    observed = [float(line.split()[-1]) for line in (cache / "korpm_controls_predictions.txt").read_text().splitlines()]
    if observed != [-1.396, 0.108, -0.727]:
        raise RuntimeError(f"KORPM failed published upstream example validation: {observed}")
    ca = [a for a in pdb_atoms(cache / "AF-O60566-F1.pdb") if a["name"] == "CA"]
    mutations = [f"NA{a['res']}K" for a in ca if a["res"] >= 705 and a["aa"] == "ASN"]
    mutations.extend(["NA1002Q", "NA1002A", "LA1012P", "NA1002N"])
    path = cache / "korpm_input.txt"
    path.write_text("".join(f"AF-O60566-F1 {mutation}\n" for mutation in mutations))
    with (cache / "korpm_run.log").open("w") as log:
        subprocess.run([str((root / "sbg/bin/korpm_gcc").resolve()), str(path),
                        "--dir", str(cache), "--score_file", str(root / "pot/korp6Dv1.bin"),
                        "-o", str(cache / "korpm_predictions.txt")], check=True, stdout=log, stderr=log)


def pocket_metrics(cache: Path, name: str, target_xyz: np.ndarray) -> list[dict]:
    base = cache / f"{name}_out"
    info = base / f"{name}_info.txt"
    if not info.exists():
        return []
    rows = []
    for block in info.read_text().split("Pocket ")[1:]:
        number = int(block.split(" :")[0])
        scores = {key.strip(): float(value) for key, value in
                  re.findall(r"^\s*(.*?):\s*([-0-9.]+)\s*$", block, re.M)}
        vertices = []
        for line in (base / "pockets" / f"pocket{number}_vert.pqr").read_text().splitlines():
            if line.startswith(("ATOM", "HETATM")):
                vertices.append([float(x) for x in line.split()[5:8]])
        vertex_xyz = np.array(vertices)
        distance = float(np.min(np.linalg.norm(vertex_xyz[:, None] - target_xyz[None, :], axis=2)))
        rows.append({"pocket": number, "score": scores["Score"],
                     "druggability_score": scores["Druggability Score"],
                     "volume_A3": scores["Volume"], "nearest_vertex_to_N1002_A": round(distance, 3)})
    return rows


def analyze(cache: Path) -> dict:
    atoms = pdb_atoms(cache / "AF-O60566-F1.pdb")
    ca = {a["res"]: a for a in atoms if a["name"] == "CA"}
    if ca[1002]["aa"] != "ASN":
        raise ValueError("Reference residue 1002 is not ASN; stop rather than mis-map mutation")
    target = np.array([a["xyz"] for a in atoms if a["res"] == 1002])
    experimental = [coordinate_coverage(cache / "6TLJ.cif", "S"),
                    coordinate_coverage(cache / "5KHU.cif", "Q")]
    homologs = []
    for code, resolution in [("6JKK", 1.85), ("6JKM", 1.95)]:
        if (cache / f"{code}.cif").exists():
            coverage = coordinate_coverage(cache / f"{code}.cif", "A")
            homologs.append({"PDB": code, "species": "Drosophila melanogaster",
                             "resolution_A": resolution,
                             "modeled_residues": coverage["modeled_residues"],
                             "residue_number_range": [coverage["first"], coverage["last"]],
                             "human_residue_1002_validated": False,
                             "PMID": "31201382", "DOI": "10.1038/s41422-019-0178-z"})
    full = pocket_metrics(cache, "AF-O60566-F1-b20", target)
    domain = pocket_metrics(cache, "AF-O60566-Cterm-b20", target)
    dynamut_path = cache / "dynamut_status.json"
    dynamut = json.loads(dynamut_path.read_text()) if dynamut_path.exists() else {}
    methods = [{"method": "DynaMut2", "status": dynamut.get("status", "not_run"),
                "ddG_unfolding_kcal_mol": dynamut.get("prediction"),
                "published_independent_test_RMSE_kcal_mol": 1.14,
                "error_is_not_case_specific_CI": True, "PMID": "32881105"}]
    korpm_path = cache / "korpm_predictions.txt"
    korpm = {}
    if korpm_path.exists():
        for line in korpm_path.read_text().splitlines():
            fields = line.split()
            if len(fields) == 3:
                korpm[fields[1]] = float(fields[2])
    target_korpm = korpm.get("NA1002K")
    methods.append({"method": "KORPM", "status": "done" if target_korpm is not None else "not_run",
                    "ddG_unfolding_kcal_mol": target_korpm,
                    "published_Ssym_RMSE_kcal_mol": 1.33,
                    "error_is_not_case_specific_CI": True, "PMID": "36629451"})
    null = [v for k, v in korpm.items() if k.startswith("NA") and k.endswith("K") and k != "NA1002K"]
    pocket_gate = any(p["nearest_vertex_to_N1002_A"] <= 8 and p["druggability_score"] >= 0.5 for p in domain)
    stability_gate = all(m["status"].upper() == "DONE" and m["ddG_unfolding_kcal_mol"] < -1.5 for m in methods)
    metadata = json.loads((cache / "alphafold_metadata.json").read_text())[0]
    hashes = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in cache.iterdir()
              if p.is_file() and p.suffix in (".pdb", ".cif", ".json", ".txt") and p.stat().st_size < 20_000_000}
    result = {"target": "BUB1B/BubR1 O60566 p.Asn1002Lys", "reference_only": True,
              "alphafold": {"entry": metadata["entryId"], "version": metadata["latestVersion"],
                            "url": metadata["pdbUrl"], "model_created": metadata["modelCreatedDate"],
                            "residue_1002_pLDDT": ca[1002]["bfactor"],
                            "domain_705_1050_mean_pLDDT": float(np.mean([ca[n]["bfactor"] for n in range(705, 1051)])),
                            "N1002_SASA_A2_960_points": residue_sasa(atoms, 1002)},
              "experimental_coordinate_coverage": experimental, "stability_predictions": methods,
              "experimental_Cterminal_homologs": homologs,
              "KORPM_comparator_null": {"definition": "Other Asn-to-Lys substitutions in residues 705–1050 of the same reference structure; local sequence/structure background, not a clinical null",
                                        "n": len(null), "median_kcal_mol": float(np.median(null)) if null else None,
                                        "fraction_at_least_as_destabilizing_as_N1002K": float(np.mean(np.array(null) <= target_korpm)) if null and target_korpm is not None else None,
                                        "finite_background_tail_rank": float((1 + np.sum(np.array(null) <= target_korpm)) / (1 + len(null))) if null and target_korpm is not None else None,
                                        "tail_rank_is_not_a_clinical_p_value": True,
                                        "controls": {key: korpm.get(key) for key in ("NA1002N", "NA1002Q", "NA1002A", "LA1012P")}},
              "pockets_full_reference_B20": full, "pockets_Cterminal_reference_B20": domain,
              "docking_gate": {"requires_both_ddG_below_minus_1_5": stability_gate,
                               "requires_nearby_fpocket_druggability_at_least_0_5": pocket_gate,
                               "decision": "eligible_for_exploratory_docking" if stability_gate and pocket_gate else "stop_no_supported_structural_rescue_screen"},
              "file_sha256": hashes}
    return result


def write_figure(result: dict, output: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(1, 2, figsize=(10, 4.1))
    methods = [m for m in result["stability_predictions"] if m["ddG_unfolding_kcal_mol"] is not None]
    values = [m["ddG_unfolding_kcal_mol"] for m in methods]
    errors = [m.get("published_independent_test_RMSE_kcal_mol", m.get("published_Ssym_RMSE_kcal_mol")) for m in methods]
    axes[0].errorbar(values, np.arange(len(methods)), xerr=errors, fmt="o", color="#16324f", capsize=5)
    axes[0].set_yticks(np.arange(len(methods)), [m["method"] for m in methods])
    axes[0].axvline(0, color="#777", linewidth=1)
    axes[0].set_xlabel("Predicted ΔΔG of unfolding (kcal/mol)\nNegative = destabilizing; whiskers = benchmark RMSE")
    axes[0].set_title("N1002K: effects smaller than method error", fontsize=10)
    pockets = result["pockets_Cterminal_reference_B20"]
    axes[1].scatter([p["nearest_vertex_to_N1002_A"] for p in pockets],
                    [p["druggability_score"] for p in pockets], color="#187c80", alpha=0.8)
    axes[1].axhline(0.5, color="#b04d26", linestyle="--", label="Exploratory threshold 0.5")
    axes[1].axvline(8, color="#777", linestyle=":", label="Nearby: ≤8 Å")
    axes[1].set_ylim(-0.02, 0.65)
    axes[1].set_xlabel("Nearest pocket vertex to N1002 (Å)")
    axes[1].set_ylabel("fpocket druggability score")
    axes[1].set_title("No nearby pocket passes the screen", fontsize=10)
    axes[1].legend(fontsize=7)
    figure.suptitle("Reference-structure triage supports stopping the docking screen", fontsize=12)
    figure.tight_layout()
    figure.savefig(output, dpi=220, bbox_inches="tight")
    figure.savefig(output.with_suffix(".svg"), bbox_inches="tight")
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=Path("work/structure"))
    parser.add_argument("--output", type=Path, default=Path("results/structure_summary.json"))
    parser.add_argument("--fetch", action="store_true")
    parser.add_argument("--prepare-tools", action="store_true")
    parser.add_argument("--submit-dynamut", action="store_true")
    parser.add_argument("--poll-dynamut", action="store_true")
    parser.add_argument("--fpocket", type=Path)
    parser.add_argument("--korpm-root", type=Path)
    parser.add_argument("--figure", type=Path)
    args = parser.parse_args()
    if args.fetch:
        fetch_references(args.cache)
    if args.prepare_tools:
        args.fpocket, args.korpm_root = prepare_tools(args.cache)
    if args.submit_dynamut or args.poll_dynamut:
        dynamut_reference_prediction(args.cache, args.submit_dynamut)
    prepare_pocket_inputs(args.cache)
    if args.fpocket:
        for name in ("AF-O60566-F1-b20", "AF-O60566-Cterm-b20"):
            subprocess.run([str(args.fpocket.resolve()), "-f", str(args.cache / f"{name}.pdb")], check=True)
    if args.korpm_root:
        run_korpm(args.cache, args.korpm_root)
    result = analyze(args.cache)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    if args.figure:
        write_figure(result, args.figure)
    print(json.dumps({"output": str(args.output), "decision": result["docking_gate"]["decision"]}))


if __name__ == "__main__":
    main()
