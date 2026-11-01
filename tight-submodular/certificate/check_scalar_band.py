#!/usr/bin/env python3
"""Independent endpoint certificate for the concise all-count manuscript band."""
from fractions import Fraction as Q
from hashlib import sha256
from pathlib import Path
import argparse
import json
import os
import sys
import tempfile

HERE = Path(__file__).resolve().parent
from flint import arb, ctx

ctx.prec = 320


def ar(value):
    q = Q(value)
    return arb(q.numerator) / q.denominator


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    checks, margins = {}, {}

    def check(name, condition):
        assert condition, name
        checks[name] = True

    def positive(name, value, threshold=0):
        margins[name] = str(value)
        check(name, value > threshold)

    z = arb(-1).exp()
    beta = 1 - z
    k = 1 / (4 * beta)
    R4 = ar("1.00032357")
    R5 = ar("1.02631852")
    K = 7 * R4.log() / 6
    bc = 43 * R4.log() / 75
    U = ar(".6905799")
    lam = ar(".2319048039")
    xi = ar(Q(1, 10**6))
    tl, th, dm = ar(".1199"), ar(".3001"), ar(Q(23, 20))
    a0 = ar(".001")
    b0 = (1 - (-10 * a0).exp()) / 10
    mu_hi = ar(Q(21, 50))
    slope = ar(Q(11, 40))
    qm, qp = R5 - th / 2, R4 - ar(Q(2, 5)) * tl
    tc = 10 * (R5 - R4)
    tf = (1 + k.log() - bc) / k

    # Basic profile, explicit credit, and root existence.
    ell = arb(2).log()
    rs = (1 / ell).log() / ell
    eps = 1 + rs - 1 / ell
    # The join-slope identity exp(ell*r_star)*ell=1 is exact by
    # the definition of r_star; only its domain inequalities need a gate.
    check("profile_logarithm_domain", 0 < ell < 1)
    check("profile_linear_piece_positive_at_join", rs > eps and eps > 0)
    check("W_small_D_is_one", th + z < 1)
    check("W_large_D_is_D", th < beta)
    check("M_at_one_identity", Q(9, 16)**2 / (1 + Q(9, 16))
          == Q(81, 400))
    check("M_at_one_uses_quadratic", Q(81, 400) < Q(5, 16))
    positive("M_HOME_branch_inactive",
             ar(Q(481, 400)).log() - ar(Q(3, 5)) * dm.log())
    check("root_exists_below_one", 1 - th > bc)
    check("root_monotone_up_to_one", th < 1)
    positive("uniform_mu_upper", 1 + mu_hi.log() - mu_hi * th - bc)
    check("floor_below_uniform_mu", k < mu_hi)
    check("root_derivative_denominator", mu_hi * th < 1)
    check("lower_floor_below_one", k < 1)

    # Price envelope: every segment has an explicit lower-bound gate.
    check("cap_order", 0 < K < a0 < ar(Q(1, 400)) < ell)
    positive("nonnegative_capped_price", b0 - K / 2)
    check("exp_quadratic_series_bound",
          Q(1, 2) / (1 - Q(1, 3000)) < Q(3, 5))
    positive("price_segment_zero_to_K",
             ar(Q(1, 2)) - slope - (5 + ar(Q(3, 5))*slope)*a0)
    positive("price_segment_K_to_a0_derivative",
             (-10*a0).exp() - slope*a0.exp())
    endpoint = ar(Q(1, 400))
    positive("price_segment_flat",
             b0 - K/2 - slope*(endpoint.exp()-1))
    for label, price in [("left", endpoint), ("right", ell)]:
        positive("price_segment_linear_" + label,
                 ar(Q(2, 5))*price-K/2-slope*(price.exp()-1))
    positive("large_price_small_Delta",
             ar(Q(2, 5))*ell-K/2-slope)

    # Physical eligibility and all-price owner bound.
    check("q0_positive", qm > ar(Q(4, 5)))
    check("q0_extrema_order", qm < qp)
    check("physical_scores_eligible", Q(23, 20) < 4*Q(4, 5))
    check("overlap_owner_sign", z > beta/2)
    check("overlap_shared_sign", beta > z)
    positive("owner_linear_slope", slope-mu_hi*beta*R5)
    check("linear_M_constant", Q(5, 16)*Q(23, 20) == Q(23, 64))
    check("large_Delta_interval_nonempty", R5 < 6-qm)
    cM = ar(Q(23, 64))

    def parallel_log(delta):
        w, y = 1 + delta/R5, qm + delta
        return w.log()+y.log()-(y+cM*w).log()

    for label, delta, lower in [
        ("left", R5, ar(".1000")),
        ("right", 6-qm, ar(".1233")),
    ]:
        positive("concave_owner_" + label,
                 parallel_log(delta)-K/2-mu_hi*beta*delta, lower)

    # A and B are simultaneous bounds; monotonicity yields finite roots.
    positive("B_h123_decreasing_mu", 1/mu_hi-th-3*beta*qp, ar(".2748"))
    check("second_score_branch_decreasing", mu_hi*R5 < ar(Q(1, 2)))
    check("crossing_location", ar(".15") < tc < th)
    cstar = ar(Q(3, 5))*ar("1.003").log()
    check("large_D_HOME_uses_flat_cap", ar("1.003").log() > a0)
    positive("large_D_HOME_localized",
             ar(Q(2, 5))*ar("1.003").log()-b0)
    check("large_D_HOME_M", ar(Q(481, 400)).log() > cstar)
    check("positive_large_D_credit", cstar > 0)

    root_brackets = [
        ("first_left", ar(".15"), arb(0), ar(".390121814249")),
        ("crossing", tc, arb(0), ar(".409249486646")),
        ("large_D_left", tl, cstar, ar(".384621435921")),
    ]
    normalized_margins = []
    for label, t, credit, root_hi in root_brackets:
        positive("root_upper_sign_" + label,
                 credit+1+root_hi.log()-root_hi*t-bc)
        check("root_upper_below_uniform_" + label, root_hi < mu_hi)
        q0 = R4-ar(Q(2, 5))*t
        for h in (1, 2, 3):
            val = -beta.log()+h*root_hi*beta*q0
            gap = U+(h-1)*lam-val
            positive(f"h{h}_normalized_{label}", gap, xi)
            normalized_margins.append(gap)

    positive("small_t_controller_below_k", 1+k.log()-k*ar(".15")-bc)
    mu_lo = k/ar("1.003")
    t = ar(".15")
    q0 = R4-ar(Q(2, 5))*t
    for h in (1, 2, 3):
        check(f"small_t_direct_increasing_h{h}",
              1-ar(Q(2*h, 5))*beta > 0)
        val = -(mu_lo*beta).log()-1+mu_lo*t+bc+h*mu_lo*beta*q0
        gap = U+(h-1)*lam-val
        positive(f"h{h}_small_t_direct", gap, xi)
        normalized_margins.append(gap)
    check("all_h123_bounded_margin_above_17e_6",
          all(x > ar(".000017") for x in normalized_margins))

    # h4 floor and active-root regions.
    check("floor_join_location", ar(".18") < tf < ar(".183") < tc)
    check("floor_segment_uses_first_score", tf < tc)
    check("HOME_linear_derivative", 1-Q(2, 5)*Q(23, 20) > 0)
    positive("HOME_M_dominates_gD",
             ar(Q(481, 400)).log()-(dm.log()-1+1/dm))
    positive("HOME_localized_quadratic",
             5-ar(Q(50, 3))*a0-ar(Q(1, 2)))
    positive("HOME_localized_join", 1-(-a0).exp()-b0)
    numer = (ar(Q(2, 5))-k)*tl-(R4-1)
    positive("h4_floor_numerator", numer)
    floor_log4_gap = numer/dm-bc
    positive("h4_floor_gap_below_log4", floor_log4_gap, ar(".0000029323763"))
    floor_gap = U+3*lam-(arb(4).log()-floor_log4_gap)
    positive("h4_floor_gap_at_rational_threshold", floor_gap, xi)
    for label, t, r in [
        ("floor_join", tf, k),
        ("crossing", tc, ar(".409249486646")),
    ]:
        q0 = R4-ar(Q(2, 5))*t
        val = -beta.log()+4*r*beta*q0
        positive("h4_active_" + label, U+3*lam-val, ar(Q(1, 50000)))

    # Direct h1 tail eliminates all external count dependencies.
    check("h1_large_head_a_max_at_zero", Q(1, 4)/Q(23, 20) > Q(1, 6))
    tail = (1+ar(Q(5, 16))*dm).log()+K/2+mu_hi*th+bc
    margins["h1_large_head_upper_bound"] = str(tail)
    check("h1_large_head_upper_below_434", tail < ar(".434"))
    positive("h1_large_head_gap", U-tail, xi)
    check("deletion_price_factor_identity",
          (1+Q(5, 16)*Q(23, 20))/6 == Q(29, 128))
    positive("large_head_deletion",
             (lam-xi).exp()-1-ar(Q(29, 128))*(K/2).exp())
    check("deletion_increment_positive", lam > xi)
    check("internal_induction_margin", Q(2, 10**6) > Q(1, 10**6))

    report = {
        "status": "PASS_CONSOLIDATED_SCALAR_BAND",
        "check_count": len(checks),
        "checks": checks,
        "directed_margins": margins,
        "domain": {"t": ["1199/10000", "3001/10000"], "D": ["0", "23/20"]},
        "counts": [1, 2, 3, 4],
        "uniform_reserve": "1/1000000",
        "precision_bits": ctx.prec,
        "positive_count_dependencies": [],
        "numerical_grid_or_optimizer": False,
        "baseline_files_read": [],
        "proof_sha256": sha256((HERE.parent/"paper"/"scalar_band.tex").read_bytes()).hexdigest(),
        "checker_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "scope": "One actual row/controller/target and prices; every 0<=b<=bcap. Direct h1 tail and internal deletion prove all counts.",
    }
    path = args.output
    with tempfile.NamedTemporaryFile("w", dir=path.parent, prefix=path.name+".", delete=False) as f:
        tmp = Path(f.name)
        json.dump(report, f, indent=2)
        f.write("\n")
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)
    print(report["status"], len(checks), "directed endpoint gates")
    print("h4 minimum:", margins["h4_floor_gap_at_rational_threshold"])
    print("h1 direct tail:", margins["h1_large_head_upper_bound"])


if __name__ == "__main__":
    main()
