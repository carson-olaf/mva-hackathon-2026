"""Freeze only public API responses used by the two documented example runs."""

import hashlib
import json
from pathlib import Path
import shutil

source = Path("work/repurpose-cache")
destination = Path("examples/public_api_cache")
destination.mkdir(parents=True, exist_ok=True)
requests = {}
for name in ("repurpose_bub1b", "repurpose_positive_control"):
    result = json.loads(Path(f"results/{name}.json").read_text())
    for item in result["api_provenance"]:
        key = hashlib.sha256(json.dumps(item["request"], sort_keys=True).encode()).hexdigest()
        requests[key] = item
for key in requests:
    shutil.copyfile(source / f"{key}.json.gz", destination / f"{key}.json.gz")
for previous in destination.glob("*.json.gz"):
    if previous.name.removesuffix(".json.gz") not in requests:
        previous.unlink()
(destination / "manifest.json").write_text(json.dumps(requests, indent=2) + "\n")
print(f"Frozen {len(requests)} public API responses; no patient variants or phenotype query payloads.")
