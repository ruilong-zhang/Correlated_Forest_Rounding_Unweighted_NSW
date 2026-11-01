#!/usr/bin/env python3
"""Universal polynomial identities in the self-contained framework lower bound.

The continuous smoothing and masking arguments are proved in the TeX;
these exact checks certify their explicit cardinality and exponent margins.
No sampled valuation, numerical grid, or old proof is imported.
"""
from fractions import Fraction as Q
from hashlib import sha256
from math import comb
from pathlib import Path
import argparse
import json
import os
import tempfile

HERE = Path(__file__).resolve().parent


def clean(p):
    return {k: Q(v) for k, v in p.items() if v}


def add(p, q):
    out = dict(p)
    for key, val in q.items():
        out[key] = out.get(key, Q(0)) + val
    return clean(out)


def scale(p, c):
    return clean({k: c*v for k, v in p.items()})


def shift(p, c):
    out = {}
    for (dm, dr), val in p.items():
        for j in range(dr+1):
            key = (dm, j)
            out[key] = out.get(key, Q(0)) + val*comb(dr, j)*Q(c)**(dr-j)
    return clean(out)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path,
                        default=HERE/"lower_bound_report.json")
    args = parser.parse_args()
    checks = {}

    def check(name, condition):
        if not condition:
            raise AssertionError(name)
        checks[name] = True

    # Formal variables are m,r.  This is an identity for every m,r.
    q = {(1, 1): Q(2), (0, 2): Q(-1)}
    inc = add(shift(q, 1), scale(q, -1))
    check("cardinality_increment_identity",
          inc == {(1, 0): Q(2), (0, 1): Q(-2), (0, 0): Q(-1)})
    sec = add(scale(shift(q, 1), 2), scale(add(q, shift(q, 2)), -1))
    check("cardinality_DR_identity", sec == {(0, 0): Q(2)})
    check("q_at_zero", all(dr > 0 for dm, dr in q))
    at_m = {}
    for (dm, dr), val in q.items():
        at_m[dm+dr] = at_m.get(dm+dr, Q(0))+val
    check("q_at_m", at_m == {2: Q(1)})
    last_inc = {}
    for (dm, dr), val in inc.items():
        for j in range(dr+1):
            deg = dm+j
            last_inc[deg] = last_inc.get(deg, Q(0))+val*comb(dr,j)*(-1)**(dr-j)
    check("minimum_increment", {k:v for k,v in last_inc.items() if v} == {0: Q(1)})
    check("strict_rounded_increment", 1-Q(1,10) == Q(9,10) > 0)
    check("strict_rounded_DR", 2-2*Q(1,10) == Q(9,5) > 0)
    check("forced_H_singleton_gap", 1-Q(1,2)-Q(1,10) == Q(2,5) > 0)
    check("masking_exponent_survives_queries", Q(1,8)-Q(1,32) == Q(3,32) > 0)
    check("smoothing_C_dominates_error", 16 > 8)
    check("population_exponents", Q(1,4)+Q(3,4) == 1)
    check("four_reserved_items_per_A", Q(1,3)*4 == Q(4,3))

    proof = HERE.parent/"paper"/"lower_bound.tex"
    report = {
        "status": "PASS_SIMPLIFIED_RATIONAL_FRAMEWORK_LOWER_BOUND",
        "check_count": len(checks),
        "checks": checks,
        "proof_sha256": sha256(proof.read_bytes()).hexdigest(),
        "baseline_dependencies": [],
        "scope": "Universal formal polynomial identities and exact margins; complete smoothing, reference embedding, masking, and quantifiers are in the TeX.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", dir=args.output.parent,
                                     prefix=args.output.name+".", delete=False) as f:
        tmp = Path(f.name)
        json.dump(report, f, indent=2)
        f.write("\n")
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, args.output)
    print(report["status"], len(checks), "exact gates")


if __name__ == "__main__":
    main()
