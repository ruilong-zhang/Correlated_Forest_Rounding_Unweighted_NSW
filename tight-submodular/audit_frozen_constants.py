#!/usr/bin/env python3
"""Check the three directed literals embedded in verify_B_simple.cpp.

These are an additional, independent exact-input/Arb check. They do not
replace the source's global coverage or analytic arithmetic arguments.
Requires python-flint 0.9.0; run without Python -O.
"""
import argparse
import json
from pathlib import Path

from flint import arb, ctx, fmpq


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    ctx.prec = 320
    denominator = 1 << 48

    def q(a, b=1):
        return arb(fmpq(a, b))

    constants = {
        "EXP_A": (281756592471775, 281756592471783, q(1, 1000).exp()),
        "CAP_B": (280072281368, 280072281376, (1-q(-1, 100).exp())/10),
        "HALF_LOG_R4": (45531063257, 45531063265,
                        q(100032357, 100000000).log()/2),
    }
    report = {"status": "PASS_FROZEN_DIRECTED_CONSTANTS",
              "precision_bits": ctx.prec, "denominator": denominator,
              "checks": {}}
    for name, (lo, hi, value) in constants.items():
        if not (q(lo, denominator) <= value <= q(hi, denominator)):
            raise AssertionError((name, value))
        report["checks"][name] = {
            "lower_integer": lo, "upper_integer": hi,
            "arb_enclosure": str(value),
        }
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(report["status"])


if __name__ == "__main__":
    main()
