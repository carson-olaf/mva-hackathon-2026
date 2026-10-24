"""Fail-closed, content-addressed completion markers for local analysis files."""
from __future__ import annotations

import hashlib
import json
import os
from functools import lru_cache
from pathlib import Path


def fingerprint(path: Path) -> dict:
    """Hash complete bytes; metadata alone cannot prove an artifact is unchanged."""
    path = Path(path)
    stat = path.stat()
    # Only an in-process optimization. Every new controller hashes files again.
    return dict(_fingerprint(str(path.resolve()), stat.st_size, stat.st_mtime_ns,
                             stat.st_ctime_ns, stat.st_ino))


@lru_cache(maxsize=256)
def _fingerprint(name, size, mtime_ns, ctime_ns, inode):
    digest = hashlib.sha256()
    path = Path(name)
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    after = path.stat()
    if (after.st_size, after.st_mtime_ns, after.st_ctime_ns, after.st_ino) != (size, mtime_ns, ctime_ns, inode):
        raise RuntimeError("Artifact changed while hashing")
    return {"sha256": digest.hexdigest(), "bytes": size}


def signature(inputs, parameters, code) -> dict:
    return {"schema": 1, "inputs": {str(p): fingerprint(Path(p)) for p in inputs},
            "parameters": parameters, "code": {str(p): fingerprint(Path(p)) for p in code}}


def atomic_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".partial")
    with temporary.open("w") as handle:
        json.dump(value, handle, indent=2)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    temporary.replace(path)


def reusable(marker: Path, expected: dict, outputs) -> bool:
    """A file, partial file or old marker alone never authorizes reuse."""
    try:
        saved = json.loads(marker.read_text())
        return (isinstance(saved, dict) and saved.get("complete") is True and saved.get("signature") == expected
                and saved.get("outputs") == {str(p): fingerprint(Path(p)) for p in outputs})
    except (OSError, ValueError, TypeError):
        return False


def complete(marker: Path, expected: dict, outputs) -> None:
    # Caller validates outputs first. The marker is the final transaction commit.
    atomic_json(marker, {"complete": True, "signature": expected,
                         "outputs": {str(p): fingerprint(Path(p)) for p in outputs}})
