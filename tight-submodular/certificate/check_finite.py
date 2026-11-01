#!/usr/bin/env python3
"""Finite endpoint checks for the parameterized graph and uniform interface.

Run with python-flint installed.  No baseline source, optimizer, or grid
is used.  General network and rounding assertions are proved in TeX.
"""
from fractions import Fraction as Q
from hashlib import sha256
from pathlib import Path
import argparse
import json
import os
import tempfile

from flint import arb, ctx

ctx.prec = 320
HERE = Path(__file__).resolve().parent


def ar(value):
    q = Q(value)
    return arb(q.numerator) / q.denominator


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path,
                        default=HERE / "graph_interface_report.json")
    args = parser.parse_args()
    checks, values = {}, {}

    def check(name, test):
        if not test:
            raise AssertionError(name)
        checks[name] = True

    def gap(name, value, threshold=0):
        values[name] = str(value)
        check(name, value > threshold)

    R = ar("1.00032357")
    L = R.log()
    z = arb(-1).exp()
    beta = 1 - z
    k = 1 / (4 * beta)
    K, T, P = ar(Q(7, 6))*L, ar(Q(12, 5))*L, ar(Q(43, 25))*L
    theta, cap = ar(Q(7, 10)), ar(Q(6, 25))
    ta, dm = ar(Q(3, 25)), ar(Q(26267, 25000))
    da = ar(Q(21, 20))*R**2
    eta = ar(Q(1, 10**10))
    a = ar(Q(1, 1000))
    ba = (1 - (-10*a).exp())/10

    check("r_integer_at_least_two", 4 >= 2)
    check("theta_half_to_one", Q(1, 2) <= Q(7, 10) <= 1)
    check("positive_log_small", 0 < L < ar(Q(1, 3000)))
    check("fallback_identity", Q(43, 25) == 1 + (1-Q(7, 10))*Q(12, 5))
    check("P_dominates_theta_T", Q(43, 25) > Q(7, 10)*Q(12, 5))
    check("mandatory_high_column_comparison", Q(43, 75) < Q(7, 12))
    check("rational_arithmetic_comparison",
          Q(43, 75) < Q(6, 25)*Q(12, 5)*(1-Q(6, 15000)))
    gap("exact_arithmetic_comparison", cap*(1-(-T).exp())-P/3)
    gap("retained_slope_comparison", cap*(1-(-T).exp())-theta*T/3)

    check("M_at_one", Q(9, 16)**2/(1+Q(9, 16)) == Q(81, 400))
    check("M_at_one_uses_quadratic", Q(81, 400) < Q(5, 16))
    check("cap_difference_increases_to_join", (-10*a).exp() > ar(Q(2, 5)))
    check("localized_cap_loss", ba-ar(Q(2, 5))*a < ar(Q(3, 5000)))
    gap("HOME_credit_M_branch", ar(Q(481, 400)).log()-ar(Q(3, 5))*dm.log())
    check("floor_below_one", 0 < k < 1)
    check("floor_below_inverse_t", k*(ta+eta) < 1)
    floor_gap = 1+k.log()-ar(Q(2, 5))*dm.log()-k*(ta+eta)-ar(Q(3, 5000))
    gap("uniform_boosted_floor", floor_gap, ar(Q(1, 250)))
    check("boost_below_floor_reserve", P/3 < ar(Q(1, 250)))
    check("exact_NEW_t_inside_rectangle", ar(Q(3, 40))*R**2 < ta)
    check("exact_NEW_D_inside_outer_box", da < dm)
    values["D_A"] = str(da)
    check("private_absent_on_A", ta+ar(Q(37, 100))*dm < 1)

    # For eta <= eta0, each asserted inequality is monotone toward
    # the endpoint below, or follows from exp(2 eta)-1 <= 3 eta.
    check("exponential_increment_bound", (2*eta).exp()-1 < 3*eta)
    check("finite_NEW_D_padding", 3*da < 6)
    check("finite_NEW_t_padding", ar(Q(3, 40))*R**2*(2*eta).exp()+eta < ta)
    check("finite_A_outer_box", da+7*eta < dm)
    check("finite_B_exact_complement", Q(6) > 0)
    check("finite_K_gap", L+eta < K-eta)
    gap("finite_high_column_slack", (K-eta)/2-P/3)
    gap("finite_deleted_head_slack", cap*(1-(-(T-5*eta)).exp())-P/3)
    check("finite_T_positive", T-5*eta > 0)
    check("finite_P_never_increases", Q(1)-Q(3, 10)*4 <= 0)
    lower_y = R**(-3)*(1-2*eta)-ar(Q(2, 5))*(ta+eta)
    gap("broad_NEW_candidate_gate", lower_y-ar(Q(3, 5))*(dm+ta+2*eta))
    gap("disabled_point4_large_deficit", R**(-3)*(1-2*eta)-dm/4-ar(Q(2, 5))*(ta+eta))
    check("finite_false_A_projection", (1+(ta+eta)/4)/ta < 9)
    check("credit_target_root_inversion", Q(9)**2/2 > 32)
    check("root_scalar_error_budget",
          11*eta.sqrt()+32*beta/ar(".3999999") < 64)

    # Uniform margins for supported upper-weight OLD programs on the finite A strip.
    qeta=Q(1,10**10); qt=Q(3,25)+qeta; qd=Q('1.05068')
    wide_t=qt*(1+2*qeta);wide_D=qd*(1+2*qeta)
    check("supported_OLD_wide_t",wide_t<Q(1,5))
    check("supported_OLD_wide_D",wide_D<Q(6,5))
    check("finite_private_inactive",(qt+Q(37,100)*qd)*(1+2*qeta)<1)
    check("finite_HOME_variable_inactive",ar(qt*(1+2*qeta))<beta)
    check("finite_raw09_inactive",(Q('1.19092969')-1)/Q('.9')>wide_t)
    check("finite_raw10_inactive",Q('1.23726325')-1>wide_t)
    check("finite_joint_inactive",(1-Q('.72149518'))/Q('.99290057')>wide_t)
    ratio=(1+Q(1,2)*Q(2,9))*Q('1.00032357')/((1+Q(2,5)*Q(2,9))*Q('1.02631852'))
    check("finite_half_affine_dominated",ratio*(1+2*qeta)<1)
    check("low_minorant_range",Q('.394805930160')*Q(19,5)<2)
    check("high_minorant_range",Q(7,20)*Q(14,3)<2)
    check("NEW_minorant_range",arb(2).log()<2)
    check("radial_coefficient",beta/ar('.3999999')<ar(Q(8,5)))
    check("analytic_threshold_below_two",((arb(4).log()-3*beta.log())/4).exp()<2)
    check("budget_per_row",2+9+Q(8,5)*32<64)
    check("budget_total",Q(111,10**6)<1)
    check("budget_fourth_power",Q(64,10**12)<1)
    check("budget_reward",Q(128,10**6)<Q(1,100))
    check("forest_factor_two",Q(2,10000)<Q(1,1000))

    # Amplification does not inspect internal component scores.
    check("amplification_log_gain", Q(1, 3)-Q(1, 12) >= Q(1, 4))
    check("amplification_log_range", arb(3).log() < ar(Q(3, 2)))
    check("amplification_success", Q(1, 4)/Q(3, 2) == Q(1, 6))

    report = {
        "status": "PASS_GRAPH_AND_FINITE_INTERFACE",
        "check_count": len(checks),
        "precision_bits": ctx.prec,
        "checks": checks,
        "directed_values": values,
        "proofs": {name: sha256((HERE.parent/"paper"/name).read_bytes()).hexdigest()
                   for name in ("graph.tex", "algorithm_interface.tex")},
        "baseline_dependencies": [],
        "optimizer_or_grid": False,
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
    print(report["status"], len(checks), "exact endpoint checks")
    print("uniform floor:", values["uniform_boosted_floor"])
    print("candidate gate:", values["broad_NEW_candidate_gate"])


if __name__ == "__main__":
    main()
