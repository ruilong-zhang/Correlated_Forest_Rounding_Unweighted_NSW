# Verification guide and manuscript map

This guide refers to the manuscript compiled on November 2, 2026. The paper's
Section 7 and Appendix F describe the current tight submodular algorithm;
Appendix E contains the subadditive numerical analysis. The programs below
verify their finite numerical obligations. They are not an end-to-end
allocation solver or a formal verification of the analytic proofs.

## Reproduce the checks

For subadditive concentration and calibration, no third-party package is needed:

```sh
cd verification
python3 verify-subadditive.py
```

Run this with Python 3.10 or later. It writes `certificate.json` and must report
`PASS: 282 concentration intervals and all constant bounds` (spacing may vary).

For the tight submodular checks, follow the dependency instructions in
[tight-submodular/README.md](tight-submodular/README.md). From that directory:

```sh
python3 check_inventory.py
python3 reproduce.py --output-dir replay-new --jobs 2
```

Use the Python environment containing `python-flint==0.9.0`, a C++17 compiler
available as `g++`, and a new output directory. Do not use Python `-O`.
`reproduce.py` runs `verify.py`, which rebuilds the elementary-function tables,
checks the numerical gates, compiles the bounded verifier, and finishes the four
covers. The wrapper then audits the frozen directed constants. The required
final statuses are `PASS_COMPLETE_SIMPLIFIED_TIGHT_PROOF_REPLAY` and
`PASS_FROZEN_DIRECTED_CONSTANTS`. A timeout, unresolved cover, assertion failure,
or nonzero exit is not a successful replay. Existing saved PASS reports are not
inputs to a fresh replay.

## Manuscript-to-code map

All code paths below are relative to the repository root. The short filenames
in the paper's Appendix F tables are relative to `tight-submodular/certificate/`.
Unless stated otherwise, saved reports in the last column are in
`tight-submodular/verification/2026-10-30/`.

| Current manuscript obligation | Verification code | Saved result |
|---|---|---|
| Appendix E.2, (53)-(54): subadditive concentration; E.3, (55): approximation constant | [verification/verify-subadditive.py](verification/verify-subadditive.py) | `verification/certificate.json` at repository root |
| F.1.1: ordinary inequality (63), joint inequality (64), private profile (66) | [local/replay_local.py](tight-submodular/certificate/local/replay_local.py) | `local.json` |
| Own-item inequality (90) and multiplier endpoint constants, Lemma F.10 | [check_home.py](tight-submodular/certificate/check_home.py) | `home.json` |
| Repair rectangle, Lemma F.13 | [check_scalar_A.py](tight-submodular/certificate/check_scalar_A.py) | `A.json` |
| Transition region, Lemma F.15 | [check_scalar_band.py](tight-submodular/certificate/check_scalar_band.py) | `band.json` |
| Bounded comparison for reserved-item counts 1-4, F.3.3 / Lemma F.16 | [verify_B_simple.cpp](tight-submodular/certificate/verify_B_simple.cpp) | `logs/h1.log` through `logs/h4.log`, summarized in `replay_report.json` |
| Unbounded singleton values, (118), and counts at least five, F.3.4 | [check_tail.py](tight-submodular/certificate/check_tail.py) | `tail.json` |
| Compatible lower bounds, F.1.3; numerical ingredients of matching repair | [check_minorants.py](tight-submodular/certificate/check_minorants.py) | `minorants.json` |
| Numerical endpoint inequalities for repair and finite implementation, F.1.4, F.2.3, F.4-F.5 | [check_finite.py](tight-submodular/certificate/check_finite.py) | `finite.json` |
| Fresh elementary-function tables and directed input data for the bounded verifier | [generate_runtime_data.py](tight-submodular/certificate/generate_runtime_data.py) | `data/runtime_data_report.json` |
| Supplementary implementation probes of bounded-verifier arithmetic and enclosures | [audit_simple_B.py](tight-submodular/certificate/audit_simple_B.py) | `primitives.json` |
| Three frozen directed constants in the C++ source | [audit_frozen_constants.py](tight-submodular/audit_frozen_constants.py) | `frozen-constants-audit.json` |

[check_menu_reductions.py](tight-submodular/certificate/check_menu_reductions.py)
checks the affine-score envelope reductions and redundancy argument used in
F.3.3 (`menu.json`). Its nondomination examples are supplementary local evidence,
not an approximation lower bound.

[check_lower_bound.py](tight-submodular/certificate/check_lower_bound.py) is an
additional inherited check. Its polynomial identities also underlie current
Lemma G.1, but its fixed margins and bound proof-source hash refer to the older
`paper/lower_bound.tex`. It does **not** certify current Appendix G wholesale;
in particular, the current rationalization uses a different rounding scale.
Its saved result is `lower_bound.json`.

## What establishes the comparison guarantee

The paper supplies the analytic reductions: it maps a fractional allocation and
its stored auxiliary optimum to an admissible scalar state, proves the local
score and repair implications, and justifies the interval coverage and tails.
The numerical routines verify the finite inequalities remaining after those
reductions. Appendix F.3.4 combines the cases to establish Lemma F.8. The
matching-repair, forest-rounding, and finite algorithm arguments then use that
lemma as specified in Section 7 and the dependency diagram.

The bounded verifier covers the two closed rectangles stated in Appendix F.3,
for each count 1, 2, 3, and 4, with the stated singleton cutoffs. It accepts a box
only using an enclosure valid throughout that box, or excludes it by a proved
infeasibility test. The transition checks and analytic tail reductions cover
the remaining cases. These checks are not a sampling experiment over scalar
states. Conversely, the supplementary arithmetic probes alone do not prove
coverage of a continuum.

`check_finite.py` checks numerical endpoint inequalities. It does not verify
arbitrary matching-repair executions or the central-point search in F.6.
Hash agreement establishes file identity, not mathematical correctness.

## Fresh upload check and historical sources

Before upload on November 2, 2026, a fresh replay completed all 4,912,892 bounded
boxes, all auxiliary gates, and the directed-constant audit. The subadditive
verifier checked all 282 intervals and reproduced its shipped certificate
byte-for-byte. [REPLAY-2026-11-02.json](REPLAY-2026-11-02.json) records the run
summary and the 33 bound source hashes. The complete October 30 reports and
logs remain in the unchanged tight-submodular release.

The 14 `.tex` files in `tight-submodular/paper/` are an immutable historical
proof snapshot, not another copy of the current paper. Some routines hash those
files for provenance; removing them without updating the provenance checks
breaks reproduction. Their older allocation-procedure description is superseded
by the current manuscript's Section 7 and Appendix F. Read this Markdown guide
and the current manuscript for the current proof interface.
