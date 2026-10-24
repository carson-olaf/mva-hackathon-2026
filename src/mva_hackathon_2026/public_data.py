"""Bounded public-reference downloads with reproducible response caching."""
import hashlib
import json
import time
from datetime import datetime, timezone
from pathlib import Path
import requests


class PublicCache:
    def __init__(self, root):
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)

    def get(self, url, *, params=None, text=False):
        key = hashlib.sha256(json.dumps([url, params], sort_keys=True).encode()).hexdigest()
        path = self.root / f"{key}.{'txt' if text else 'json'}"
        if path.exists():
            return path.read_text() if text else json.loads(path.read_text())
        error = None
        for attempt in range(4):
            try:
                response = requests.get(url, params=params, headers={"Accept": "application/json", "Content-Type": "application/json"}, timeout=90, stream=True)
                response.raise_for_status()
                if int(response.headers.get("content-length", 0)) > 1_000_000_000:
                    raise ValueError("Downloads over 1 GB require user approval")
                chunks, size = [], 0
                for chunk in response.iter_content(1024 * 1024):
                    size += len(chunk)
                    if size > 1_000_000_000:
                        raise ValueError("Download exceeded approved size limit")
                    chunks.append(chunk)
                raw = b"".join(chunks)
                value = raw.decode() if text else json.loads(raw)
                path.write_text(value if text else json.dumps(value))
                path.with_suffix(".meta.json").write_text(json.dumps({"url": response.url, "retrieved_utc": datetime.now(timezone.utc).isoformat(), "sha256": hashlib.sha256(raw).hexdigest(), "bytes": size}, indent=2))
                return value
            except requests.RequestException as exc:
                error = exc
                time.sleep(2 ** attempt)
        raise error
