#!/usr/bin/env python3
"""Check distributed file integrity; this does not replay numerical proofs."""
from hashlib import sha256
import json
from pathlib import Path

def main():
    root = Path(__file__).resolve().parent
    inventory = json.loads((root / "release-inventory.json").read_text())
    for name, record in inventory["files"].items():
        data = (root / name).read_bytes()
        if len(data) != record["bytes"] or sha256(data).hexdigest() != record["sha256"]:
            raise SystemExit("FAIL_FILE_INTEGRITY " + name)
    sources = json.loads((root / "provenance/verified-source-sha256.json").read_text())
    canonical = json.loads((root / "provenance/canonical-manifest.json").read_text())
    for name, expected in sources.items():
        if sha256((root / name).read_bytes()).hexdigest() != expected:
            raise SystemExit("FAIL_SOURCE_HASH " + name)
        if canonical["files"][name]["sha256"] != expected:
            raise SystemExit("FAIL_CANONICAL_SOURCE_HASH " + name)
    exports = json.loads((root / "provenance/report-export.json").read_text())
    for name, record in exports["files"].items():
        if sha256((root / name).read_bytes()).hexdigest() != record["distributed_sha256"]:
            raise SystemExit("FAIL_REPORT_EXPORT_HASH " + name)
    print("PASS_RELEASE_INVENTORY", len(inventory["files"]), "files;", len(sources), "unchanged source inputs")

if __name__ == "__main__":
    main()
