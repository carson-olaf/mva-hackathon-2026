"""Small public API requests with an auditable, content-checked offline cache."""

from datetime import datetime, timezone
import gzip
import hashlib
import json
import os
from pathlib import Path
import tempfile
import time

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


class PublicCache:
    def __init__(self, root: Path | str, offline: bool = False, fallbacks: tuple[Path, ...] = ()):
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)
        self.offline = offline
        self.fallbacks = fallbacks
        self.used: dict[str, dict] = {}
        self.session = requests.Session()
        self.session.headers["User-Agent"] = "mva-repurpose/0.1 public-gene-drug-research"
        self.session.mount("https://", HTTPAdapter(max_retries=Retry(
            total=3, backoff_factor=0.8, status_forcelist=[429, 500, 502, 503, 504],
            allowed_methods=["GET", "POST"], respect_retry_after_header=True,
        )))

    def json(self, url: str, *, params: dict | None = None, body: dict | None = None):
        request = {"url": url, "params": params or {}, "body": body}
        key = hashlib.sha256(json.dumps(request, sort_keys=True).encode()).hexdigest()
        path = self.root / f"{key}.json.gz"
        if not path.exists():
            path = next((p / path.name for p in self.fallbacks if (p / path.name).exists()), path)
        if path.exists():
            with gzip.open(path, "rt") as handle:
                record = json.load(handle)
        else:
            if self.offline:
                raise FileNotFoundError(f"Offline public cache miss: {url}; key={key}")
            time.sleep(0.15)
            response = self.session.request(
                "POST" if body is not None else "GET", url, params=params,
                json=body, timeout=(15, 120),
            )
            response.raise_for_status()
            if len(response.content) > 25_000_000:
                raise ValueError("Public API response exceeds the 25 MB safety limit")
            payload = response.json()
            if isinstance(payload, dict) and payload.get("errors"):
                raise ValueError(f"Public API returned errors: {payload['errors']}")
            record = {
                "request": request, "fetched_at": datetime.now(timezone.utc).isoformat(),
                "payload": payload, "http_status": response.status_code,
            }
            record["payload_sha256"] = self._digest(payload)
            descriptor, temporary = tempfile.mkstemp(dir=self.root, suffix=".tmp")
            os.close(descriptor)
            try:
                with gzip.open(temporary, "wt") as handle:
                    json.dump(record, handle, sort_keys=True)
                os.replace(temporary, path)
            finally:
                Path(temporary).unlink(missing_ok=True)
        if record["request"] != request or self._digest(record["payload"]) != record["payload_sha256"]:
            raise ValueError(f"Public API cache integrity check failed: {path}")
        self.used[key] = {k: v for k, v in record.items() if k != "payload"}
        return record["payload"]

    @staticmethod
    def _digest(payload) -> str:
        return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()

    def pages(self, resource: str, collection: str, params: dict):
        url = f"https://www.ebi.ac.uk/chembl/api/data/{resource}.json"
        current = self.json(url, params={**params, "limit": 1000})
        values = list(current[collection])
        while current["page_meta"].get("next"):
            next_url = current["page_meta"]["next"]
            if next_url.startswith("/"):
                next_url = "https://www.ebi.ac.uk" + next_url
            current = self.json(next_url)
            values.extend(current[collection])
        return values
