#!/usr/bin/env python3
"""Independent exact constants for the consolidated repair-rectangle proof.

Only rational arithmetic and directed Arb elementary functions are used.
There is no parameter grid, optimization routine, or historical proof import.
"""
from fractions import Fraction as Q
from hashlib import sha256
from pathlib import Path
import argparse
import json
import os
import sys
import tempfile

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
from flint import arb, fmpq, ctx
ctx.prec = 320


def ar(x):
    x = Q(x)
    return arb(fmpq(x.numerator, x.denominator))


def M(x):
    x = Q(x)
    v = Q(9, 16) * x
    return 1 + min(v * v / (1 + v), Q(5, 16) * x)


def atomic_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent,
                                     prefix=path.name + ".", suffix=".tmp",
                                     delete=False) as out:
        json.dump(value, out, indent=2)
        out.write("\n")
        out.flush()
        os.fsync(out.fileno())
        tmp = out.name
    os.replace(tmp, path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=HERE / "scalar_A_report.json")
    args = parser.parse_args()
    checks = []

    def pos(name, value):
        assert value > 0, (name, str(value))
        checks.append({"name": name, "strict_margin": str(value)})

    z = arb(-1).exp()
    beta = 1 - z
    k = 1 / (4 * beta)
    R0, R4 = Q("1.00024010"), Q("1.00032357")
    T, Dm, t0, D0 = Q(3, 25), Q(26267, 25000), Q(9, 125), Q(101, 100)
    a0, eps, reserve = Q(1, 1000), Q(3, 5000), Q(1, 250)
    L = ar(R4).log()
    K, bc = 7 * L / 6, 43 * L / 75
    zz = 1 - (-12 * L / 5).exp()
    b0 = (1 - (-10 * ar(a0)).exp()) / 10
    U = (arb(4).log() - 3 * beta.log()) / 4
    lam = (4 * beta).log() / 4
    vartheta = k / ar(Q(2, 5))

    def gd(x):
        return ar(x).log() - 1 + 1 / ar(x)

    # The one floor controller and the ordinary owner envelope.
    pos("canonical_A_contained", Dm - Q(21, 20) * R4 ** 2)
    pos("unconditional_slope_above_k", ar("0.395494177") - k)
    pos("floor_admissible", 1 - k * ar(T))
    pos("local_cap_loss_maximum_at_a0", (-10 * ar(a0)).exp() - ar(Q(2, 5)))
    pos("local_HOME_loss_below_eps", ar(eps) - b0 + ar(Q(2, 5) * a0))
    assert M(1) == Q(481, 400)
    pos("M1_above_Dmax", M(1) - Dm)
    pos("cap_inside_smooth_piece", ar(a0) - K)
    pos("plateau_price_pays_half_cap", b0 - K / 2)
    pos("small_price_nonnegative", Q(1, 2) - 5 * a0)
    floor_H = 1 + k.log() - ar(Q(2, 5)) * ar(Dm).log() - k * ar(T) - ar(eps)
    pos("h0_floor_reserve", floor_H - ar(reserve))
    pos("floor_exceeds_any_boost", ar(reserve) - bc)
    pos("h0_missing_heads_paid", ar(reserve) - 4 * bc)
    pos("shared_gamma_coefficient_for_k", z - ar(Q(1, 4)))
    pos("shared_gamma_coefficient_for_point4", z - beta * ar(Q(2, 5)))
    pos("owner_positive_denominator", ar(R0) - k * ar(T))
    pos("owner_argument_below_30_over_23", ar(Q(30, 23)) * (ar(R0) - k * ar(T)) - ar(Dm))
    pos("owner_M_below_four_thirds", Q(4, 3) - M(Q(30, 23)))
    pos("power_chord_at_two", ar(Q(2, 5)) * arb(2).log() - ar(R0 / 4))
    pos("linear_chord_at_two", ar(Q(3, 2)).log() - ar(R0 / 4))
    pos("linear_chord_at_eight", arb(6).log() - ar(7 * R0 / 4))

    count_bounds = {}
    for h in (1, 2, 3):
        ch = ar(Q(h, 4) * R0) + ar(1 - Q(h, 4)) * k * ar(T)
        if h == 1:
            pos("h1_D_derivative_positive", ar(Q(2, 5)) - ch)
            upper = (arb(4).log() - 1 + ar(Q(2, 5)) * ar(Dm).log()
                     + ch / ar(Dm) + ar(eps) + K / 2)
        else:
            pos(f"h{h}_D_derivative_negative", ch - ar(Q(2, 5) * Dm))
            upper = arb(4).log() - 1 + ch + ar(eps) + h * K / 2
        margin = U + (h - 1) * lam - upper
        pos(f"h{h}_bounded_reserve", margin - ar(reserve))
        charge = (4 - h) * bc + ar(Q(6 * h, 25)) * zz
        pos(f"h{h}_all_extra_charges_paid", ar(reserve) - charge)
        count_bounds[str(h)] = {"bounded_upper": str(upper),
                                "bounded_reserve": str(margin),
                                "extra_charge_upper": str(charge)}
    pos("h1_tail_max_at_A_zero", (8 * ar(R0) - k * ar(T)) / (4 * ar(Dm)) - 1)
    tail1 = ar(Q(4, 3) * R0).log() + k * ar(T) + ar(eps) + K / 2
    pos("h1_unbounded_tail_reserve", U - tail1 - ar(reserve))
    inv_tail = (K / 2).exp() / 6
    pos("tail_inverse_price_below_0167", ar(Q(167, 1000)) - inv_tail)
    pos("lambda_exp_above_five_fourths", 4 * beta - ar(Q(5, 4) ** 4))
    pos("deletion_rational_margin", Q(5, 4) - 1 - Q(167, 1000))
    for h in (2, 3):
        pos(f"h{h}_same_state_deletion", lam - (1 + inv_tail / (h - 1)).log())

    # The tight, NEW-supported low-price edge.
    cut = Q(199, 200)
    pos("critical_convex_weight_below_one", 1 - vartheta)
    pos("critical_NEW_mean_eligibility", ar(Q(3, 40)) * (ar(cut) - k * ar(t0)) - ar(t0))
    pos("critical_NEW_scale_eligibility", ar(Q(21, 20)) * (ar(cut) - k * ar(t0)) - ar(D0))
    pos("critical_tiny_head_pays_arithmetic", ar(Q(1, 200)) - ar(4 * D0 * Q(6, 25)) * zz)
    pos("critical_exponential_quadratic", ar(Q(6, 5)) - (2 * ar(a0)).exp())
    pos("critical_lowprice_slope", 1 - vartheta - ar(Q(112, 5)) * K)
    pos("critical_arithmetic_coefficient", vartheta - ar(Q(24, 25) * D0))
    pos("critical_middle_derivative", 4 * (-10 * ar(a0)).exp() - 2 * (2 * ar(a0)).exp())
    pos("critical_plateau_endpoint", 1 + 4 * b0 - 2 * K - ar(R0) * ar(Q(1, 400)).exp())
    for point, name in ((ar(Q(1, 400)), "left"), (arb(2).log(), "right")):
        pos("critical_concave_price_" + name, 1 + ar(Q(8, 5)) * point - 2 * K - ar(R0) * point.exp())

    # The conditional ordinary edge on the remaining time strip.
    pos("strip_theta_above_arithmetic_slope", Q(1, 4) / D0 - Q(6, 25))
    pos("strip_tiny_head_paid", ar((R4 - Q(9, 10)) / (4 * D0)) - ar(Q(6, 25)) * zz)
    pos("strip_ordinary_eligibility", 4 * (Q(9, 10) - Q(2, 5) * T) - Dm)
    pos("strip_lowprice_derivative", (-10 * K).exp() - ar(Q(1, 2)) - ar(R4) * K.exp() / 4 - ar(Q(6, 25)))
    pos("strip_middle_derivative", (-10 * ar(a0)).exp() - ar(R4) * ar(a0).exp() / 4)
    exp_edge = ar(R4) * (ar(Q(1, 400)).exp() - 1) / 4
    pos("strip_plateau_endpoint", b0 - K / 2 - exp_edge)
    pos("strip_linear_left", ar(Q(2, 5 * 400)) - K / 2 - exp_edge)
    pos("strip_linear_right", ar(Q(2, 5)) * arb(2).log() - K / 2 - ar(R4 / 4))
    time_margin = (ar(Q(2, 5)) - k) * ar(t0) - ar(R4 - 1)
    pos("exact_0072_time_seam", time_margin)

    # One deletion from the proved A3 ledger replaces the historical C3 proof.
    highprice_margin = (3 * lam.exp() - 3 * (ar(Q(6, 25)) * zz).exp()
                        - (-ar(Q(2, 5)) * arb(2).log() + K / 2).exp())
    pos("A4_highprice_deleted_from_A3", highprice_margin)

    # One unconditional bridge covers all t<=.12 above D=1.01.
    pos("global_chord_endpoint", Q(8, 5) ** 2 - Q(6, 5) ** 5)
    pos("global_product_chord_coefficient", 4 * (ar(R0) - k * ar(T)) - 3 * ar(R0))
    pos("HOME_linear_difference_increases", 1 - Q(2, 5) * Dm)
    pos("HOME_bridge_above_local_onset", ar(D0).log() - ar(a0))
    pos("HOME_M_base_above_g", ar(M(1)).log() - gd(Dm))
    Cstar = 2 * K + ar(Q(24, 25)) * zz + bc
    budget = ar(R0 - 1) + Cstar
    bridges = {
        "multiplicative": ar(M(1)).log() - gd(Dm) - budget,
        "linear": ar(Q(3, 5)) * ar(D0).log() - gd(D0) - budget,
        "localized": 1 - 1 / ar(D0) - b0 - budget,
    }
    for name, margin in bridges.items():
        pos("full_upper_scale_bridge_" + name, margin)

    proofs = [ROOT / "paper/scalar_A.tex", ROOT / "paper/controller.tex",
              ROOT / "paper/scalar_state.tex"]
    report = {
        "status": "PASS_CONSOLIDATED_SCALAR_A", "exact_gate_count": len(checks),
        "precision_bits": ctx.prec, "parameter_grid": False,
        "historical_proof_imports": [],
        "scope": {"t": ["0", str(T)], "D": ["0", str(Dm)],
                  "counts": [0, 1, 2, 3, 4], "prices": "all nonnegative",
                  "controller": "actual uniform boosted controller, proved equal to k/max(1,D)"},
        "payment": "half capped prices + clipped low-price arithmetic + (4-h)bc; extra bc for h4,D>=1.01",
        "floor_H_lower": str(floor_H), "count_bounds": count_bounds,
        "h1_tail_upper": str(tail1), "exact_time_seam_margin": str(time_margin),
        "highprice_deletion_margin": str(highprice_margin),
        "bridge_margins": {key: str(value) for key, value in bridges.items()},
        "checks": checks,
        "sha256": {str(path.relative_to(ROOT)): sha256(path.read_bytes()).hexdigest()
                   for path in proofs + [Path(__file__)]},
    }
    atomic_json(args.output, report)
    print(report["status"], len(checks), "exact gates")
    print("floor H:", floor_H)
    print("time seam:", time_margin)
    print("A4 highprice deletion:", highprice_margin)
    print("bridge:", json.dumps(report["bridge_margins"]))


if __name__ == "__main__":
    main()
