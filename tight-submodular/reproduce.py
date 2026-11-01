#!/usr/bin/env python3
"""Run the unchanged numerical replay and additional frozen-constant audit."""
import argparse
from pathlib import Path
import subprocess
import sys

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--jobs", type=int, choices=(1, 2), default=2)
    args = parser.parse_args()
    if sys.flags.optimize:
        parser.error("Run Python without -O")
    root = Path(__file__).resolve().parent
    output = args.output_dir.resolve()
    if output.exists():
        parser.error("Use a new output directory")
    subprocess.run([sys.executable, str(root / "verify.py"),
                    "--output-dir", str(output), "--jobs", str(args.jobs)], check=True)
    subprocess.run([sys.executable, str(root / "audit_frozen_constants.py"),
                    "--output", str(output / "frozen-constants-audit.json")], check=True)

if __name__ == "__main__":
    main()
