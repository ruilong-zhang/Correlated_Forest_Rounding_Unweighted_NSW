SIMPLIFIED TIGHT MATCHING--REMATCHING PROOF PACKAGE
19 October 2026

READING AND SCOPE

Start with nsw_tight_simplified.pdf. The English manuscript is self-contained:
it states the algorithm, all mathematical reductions, the common rounding,
the finite rational value-oracle implementation, and the framework lower bound.
It uses the same 10-point type, A4 page size and 22 mm margins as the frozen
86-page predecessor; the simplified manuscript is 34 pages.

The exact infimum factor is

    R_star = (4 / (1 - exp(-1))^3)^(1/4)
           = 1.9948721225606033369...

For every positive rational epsilon, the finite algorithm achieves
R_star + epsilon in expected Nash social welfare. Repetition with selection
by actual welfare gives the corresponding high-probability guarantee.
This is tightness within the reserved-item matching--allocation--rematching
framework. It does not claim a finite zero-tolerance endpoint or a universal
information-theoretic lower bound for all submodular NSW algorithms.

WHAT WAS SIMPLIFIED

Five affine labels become two piecewise-affine physical programs, joined by
one joint program. This groups programs with the same preserved function;
the five affine pieces and eight protected mixture factors are retained.
The physical copies of private guards are redundant and removed.

Every row uses one boosted HOME controller. Only the selected residual
iterate is classified by one rectangle. The final algorithm samples one
matching, solves one minimum-cost flow, and uses one common forest rounding.
It does not repair every matching component or compare unknown repaired means.

One pointwise completion theorem replaces the old hybrid charge cases. Its
unit coefficient follows from rational inequalities. A single endpoint proof
serves the controller. The A proof derives four-head high-price deletion from
its already-paid three-head case and uses one upper-D bridge. An independent
analytic band treats all counts one through four. One directed C++ verifier
handles the remaining B cases, without NEW square-root or physical-guard terms.

RECOMPUTE ALL COMPUTATIONAL OBLIGATIONS

Dependencies: Python 3 with python-flint, and g++ with C++17 support. The
published replay used Python 3.12.14 and python-flint 0.9.0 on x86-64 Linux.
The exact compiler version and invocation are in verification/build_receipt.json.
No third-party library, executable, generated transcendental table, old source
package, optimizer, or previously generated cover is bundled as a prerequisite.

From the unpacked package directory, with those dependencies installed:

    python3 verify.py --output-dir fresh_replay --jobs 4

Use a NEW output directory each time. Jobs may be 1, 2, 3 or 4. Do not use
Python -O. The driver compiles with assertions enabled and generates all
runtime tables directly with Arb. The fresh complete run takes a few minutes
on the recorded environment; runtime varies with hardware and concurrency.

Success requires the exact terminal status

    PASS_COMPLETE_SIMPLIFIED_TIGHT_PROOF_REPLAY

The driver first reruns all local and analytic gates, then compiles and audits
the actual C++ arithmetic, and launches all four bounded covers. A cover must
exit successfully, exhaust both roots, have one COMPLETE terminal, satisfy the
binary-tree counting identity and exact domain/threshold checks, and leave its
source, generated tables, rational inputs and executable unchanged. The final
gate also checks that every mathematical source and certificate input is
unchanged from launch. The supplied reports are provenance; the command never
uses an existing PASS report as a proof prerequisite.

CERTIFICATE MAP

paper/                  Complete LaTeX manuscript source.
verify.py               Single entry point for the complete numerical replay.
certificate/local/      Six files: exact protected-mixture parameters, one
                        private profile, arithmetic and geometric reductions,
                        and their fresh certificate replay.
certificate/check_*     Completion/minorants, HOME, A, analytic band, finite
                        graph interface, lower bound and tail constant gates.
certificate/scalar.py   Compact shared Arb HOME and scalar definitions.
certificate/generate_runtime_data.py
                        Rebuilds all directed transcendental nodes and inputs.
certificate/verify_B_simple.cpp
                        One 48-bit directed whole-box verifier for B1--B4.
certificate/audit_simple_B.py
                        Independent arithmetic/credit/root implementation audit.
verification/           The new complete replay, current reports, launch/build
                        receipts and logs; no historical report chains.
manifest.json           SHA-256 and size of every other shipped file.

The local replay covers all 4,001 private profile nodes, 727 unconditional
ray checks, 37,304 conditional checks and 64 joint planar vertices. The compact
HOME proof covers all 4,000 mean cells with 13,806 accepted q intervals. The
remaining bounded B covers visit 4,912,892 boxes: 1,078,272 for h=1; 2,061,330
for h=2; 1,441,338 for h=3; and 331,952 for h=4. Large coordinates, zero heads,
and all larger counts are handled by the stated analytic and directed-tail
reductions. The private interpolation, rounding directions, actual-state
coupling and freshly executed report handoffs were cross-checked separately.

These are mathematical proofs assisted by exact rational and interval
computation. General lemmas and oracle complexity are proved in the manuscript;
the replay is not a proof-assistant formalization. The primitive probe audit
is an implementation check; sampled probes do not replace whole-box coverage.

REBUILD THE PDF (OPTIONAL)

With a LaTeX installation containing the packages in paper/main.tex, run
pdflatex on main.tex from paper/ three times. The compiled PDF is already
included. The certificate replay does not depend on LaTeX.
