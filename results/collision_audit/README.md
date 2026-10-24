# Public target-identity audit: reproduction

This directory contains public-derived audit findings only. Raw public inputs stay in ignored `work/public-collision-audit/`. No patient or gated inputs are needed.

Main result: the DGIdb 2024-12 export contains **13,146 wrong gene–GtoP-ligand pairs**, or **12,994 gene–normalized-drug pairs**, across **1,403 genes**. The exact numeric-swap count in the tested 2026-10 GtoP subset is zero. Eleven different current identity conflicts remain. Read the [report](../../report/collision_audit.md) for definitions and limits.

From a checkout or the extracted public replay archive:

```sh
python3 scripts/replay_collision_audit.py
python3 scripts/independent_collision_check.py
```

The first command is offline, uses only the standard library and writes to `work/collision-replay-results/`. It verifies 52 pinned inputs and eight byte-identical outputs. The independent checker uses a separate implementation and reports the main December/June/October counts and Orphanet exposure join. Python 3.12 was used for validation.

If public inputs are missing, install the project's dependencies (`uv sync`), then explicitly request public downloads:

```sh
uv run python scripts/replay_collision_audit.py --download
```

This requires `requests`. The bounded fetcher caps received response bodies at 250 MB and records request/response metadata. Replay inputs total **140,771,610 bytes**. Fresh retrieval receipts remain separate from the original acquisition metadata in `input_manifest.json`. Some source URLs are mutable: if their bytes change, reproduction **fails**. Use the preserved public archive for exact historical replay; do not erase a mismatch and call a fresh run the original result.

Behavioral tests use only synthetic fixtures:

```sh
UV_CACHE_DIR=work/.uv-audit-cache uv run --offline pytest tests/test_collision_audit.py
```

Files:

- `summary.json`: denominators, categories, source versions, disease exposure and input hashes.
- `identity_conflicts.csv`: observed conflicts across four snapshots. Filter snapshot `2024-12` and classification `confirmed_numeric_namespace_collision` for the headline count; do not add repeated snapshots together.
- `current_identity_conflicts.csv`: 11 current conflicts of an unestablished cause.
- `named_cases.csv`: seven representative cases, exact IDs and opened pharmacology references.
- `rare_disease_exposure.csv`: all 512 ORPHA entities, including 140 subtypes. The ranking measures erroneous candidate leads, not patient harm.
- `rare_disease_gene_evidence.csv`: original Orphanet association validation identifiers. These source-supplied references were not all individually opened.
- `engine_source_audit.json`: complete ChEMBL mechanism scan and bounded Open Targets identity checks.
- `live_api_evidence.json`: current DGIdb sentinel evidence and request provenance.
- `input_manifest.json` and `expected_outputs.json`: pinned inputs and exact expected result hashes.

Approval fields are database flags, not independent regulator verification. `source_pmids` preserves references supplied by GtoP; only `opened_pharmacology_pmid` in the named-case file indicates a reference inspected for this audit.

The raw cache is not tracked in Git. `python3 scripts/package_collision_audit.py` creates a local public-only replay ZIP from an explicit allowlist. It does not publish or upload anything. Source-specific terms are in [NOTICE.md](NOTICE.md).

The final [validation receipt](validation.json) records an extracted-archive replay in a fresh temporary directory with network connections disabled, the test result and the archive hash. The receipt stays outside the archive to avoid a self-referential hash.
