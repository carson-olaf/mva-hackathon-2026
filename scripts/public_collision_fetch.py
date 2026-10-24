"""Bounded, content-addressed cache for the public namespace audit.

No patient files are accepted. Network activity is explicit; analysis uses saved
responses. The 250 MB ceiling is the user's approval limit for this audit.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from urllib.parse import urlsplit

import requests

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "work/public-collision-audit"
LIMIT = 250_000_000
HOSTS = {
    "dgidb.org", "www.dgidb.org", "beta.dgidb.org", "api.dgidb.org",
    "api.github.com", "raw.githubusercontent.com", "github.com", "media.githubusercontent.com",
    "www.guidetopharmacology.org", "www.guidetoimmunopharmacology.org",
    "www.guidetomalariapharmacology.org", "www.orphadata.org",
    "www.orphadata.com", "sciences.orphadata.com", "ftp.ncbi.nlm.nih.gov",
    "ftp.ebi.ac.uk", "www.ebi.ac.uk", "api.platform.opentargets.org",
    "reactome.org", "rest.ensembl.org", "storage.googleapis.com",
    "ftp.expasy.org", "api.orphadata.com", "rest.genenames.org",
}


def fetch(url: str, name: str, *, body=None, cap=80_000_000):
    if urlsplit(url).hostname not in HOSTS or not url.startswith("https://"):
        raise ValueError("URL must use HTTPS and an approved public host")
    if Path(name).name != name or name in {"manifest.json", "state.md"}:
        raise ValueError("Use a plain cache filename")
    CACHE.mkdir(parents=True, exist_ok=True)
    manifest_path = CACHE / "manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {
        "approved_limit_bytes": LIMIT, "requests": [],
    }
    path = CACHE / name
    existing = [r for r in manifest["requests"] if r.get("filename") == name and r.get("complete")]
    if path.exists() and existing:
        record = existing[-1]
        if record["url"] != url or record.get("request_body") != body:
            raise ValueError("Cache filename already belongs to another request")
        if hashlib.sha256(path.read_bytes()).hexdigest() != record["sha256"]:
            raise ValueError("Cache hash mismatch")
        return record
    used = sum(r.get("received_bytes", 0) for r in manifest["requests"])
    available = min(cap, LIMIT - used)
    if available <= 0:
        raise ValueError("Approved byte budget exhausted")
    record = {"url": url, "request_body": body, "filename": name,
              "retrieved_utc": datetime.now(timezone.utc).isoformat(),
              "complete": False, "received_bytes": 0}
    part = path.with_suffix(path.suffix + ".part")
    try:
        with requests.request("POST" if body is not None else "GET", url,
                json=body, stream=True, timeout=(20, 60),
                headers={"User-Agent": "public-target-identity-audit/1.0", "Accept-Encoding": "identity"}) as response:
            record.update(status=response.status_code, final_url=response.url,
                          content_type=response.headers.get("Content-Type"),
                          last_modified=response.headers.get("Last-Modified"),
                          etag=response.headers.get("ETag"))
            response.raise_for_status()
            if "login" in urlsplit(response.url).path.lower():
                raise ValueError("Login required; will not cache a login page as data")
            if int(response.headers.get("Content-Length", 0)) > available:
                raise ValueError("Response exceeds remaining approved byte budget")
            digest = hashlib.sha256()
            with part.open("wb") as handle:
                for chunk in response.iter_content(65536):
                    record["received_bytes"] += len(chunk)
                    if record["received_bytes"] > available:
                        raise ValueError("Response exceeded byte cap")
                    digest.update(chunk)
                    handle.write(chunk)
            part.replace(path)
            record.update(complete=True, sha256=digest.hexdigest())
    except Exception as error:
        part.unlink(missing_ok=True)
        record["error"] = f"{type(error).__name__}: {error}"
        raise
    finally:
        manifest["requests"].append(record)
        temp = manifest_path.with_suffix(".tmp")
        temp.write_text(json.dumps(manifest, indent=2) + "\n")
        temp.replace(manifest_path)
    return record


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("url")
    parser.add_argument("name")
    parser.add_argument("--body", help="JSON request body")
    parser.add_argument("--cap", type=int, default=80_000_000)
    args = parser.parse_args()
    print(json.dumps(fetch(args.url, args.name, body=json.loads(args.body) if args.body else None, cap=args.cap)))
