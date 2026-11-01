# NSW forest-rounding numerical certificates

This repository verifies the numerical proof obligations used in the tight
submodular Nash social welfare analysis. It is **not an end-to-end NSW allocation
solver**. The allocation algorithm, structural lemmas, rounding guarantees, and
the analytic reduction to these obligations are proved separately in the paper.

The verifier uses exact rational inputs, directed interval arithmetic, fresh Arb
tables, and exhaustive subdivision. Numerical PASS results do not replace the
analytic proofs that justify the domains, local score inequalities, matching
repair, or oracle implementation. In particular, the gate named
`PASS_GRAPH_AND_FINITE_INTERFACE` checks numerical endpoint obligations, not
arbitrary flow networks or full algorithm executions.

## Reproduce

Requirements: Python 3.12, `python-flint==0.9.0`, and a C++17 compiler available as
`g++`. The recorded run used Python 3.12.14, python-flint 0.9.0, and Apple clang
21.0.0 on arm64 macOS. Other platforms have not been freshly tested for this
release. The verifier compiles its C++ source with assertions enabled. No
precompiled executable, virtual environment, or generated runtime table is
distributed.

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python check_inventory.py
.venv/bin/python reproduce.py --output-dir replay-new --jobs 2
```

Use a new output directory for each run. Do not use Python `-O`. The wrapper
first calls the authoritative `verify.py`, then performs the additional frozen
constant audit. `--jobs 1` reduces local concurrency. The complete recorded run
took about three minutes; runtime depends on hardware.

The underlying command, without the additional frozen constant audit, is:

```sh
.venv/bin/python verify.py --output-dir replay-new --jobs 2
```

The complete driver must terminate with
`PASS_COMPLETE_SIMPLIFIED_TIGHT_PROOF_REPLAY`. The additional audit must terminate
with `PASS_FROZEN_DIRECTED_CONSTANTS`. Existing stored reports are never used as
prerequisites for a new proof replay.

## Recorded verification

The fresh replay completed on 30 October 2026 (Asia/Shanghai). All gates passed,
all 33 mathematical/certificate source hashes were unchanged, and all four
bounded covers exhausted their roots:

| Reserved-item count | Visited boxes | Accepted leaves | Excluded leaves |
|---:|---:|---:|---:|
| 1 | 1,078,272 | 335,271 | 203,866 |
| 2 | 2,061,330 | 813,856 | 216,810 |
| 3 | 1,441,338 | 597,214 | 123,456 |
| 4 | 331,952 | 145,772 | 20,205 |
| **Total** | **4,912,892** | **1,892,113** | **564,337** |

The run also checked all 4,001 private-profile nodes, 727 unconditional rays,
37,304 conditional checks, 64 joint planar vertices, 4,000 HOME intervals,
the analytic-region constant gates, and a 411-box singleton-tail cover. The
arithmetic/credit/root probe audit supplements the soundness review; finite
probes alone do not certify a continuum. Large coordinates, zero heads, and
larger integer head counts require the analytic reductions in the paper.

The included report copies and logs are in `verification/2026-10-30/`. Only
personal absolute path prefixes have been normalized. Their numerical values
and source hashes are unchanged. Embedded hashes refer to the original replay
artifacts; `provenance/report-export.json` gives both original and distributed
hashes. The original complete replay report SHA-256 is
`69146327c557b8390a3ad2c85d4b993ed21c377510ab0a3c569928ff0a2454be`.

## Files and provenance

- `verify.py`: authoritative, unchanged numerical replay driver.
- `certificate/`: complete verifier source and exact rational certificate
  inputs, including the private profile. Fresh transcendental tables are
  generated during every run.
- `audit_frozen_constants.py`: additional 320-bit Arb check of three hard-coded
  directed constants. The constants and their correspondence to the C++ source
  were independently inspected.
- `paper/`: frozen inherited proof-source snapshot, retained byte-for-byte
  because its hashes are part of the replay provenance. **This snapshot is not
  the current manuscript or the current central-cut algorithm description.**
  In particular, its older finite-procedure discussion should not be read as
  the current allocation algorithm. Consult Section 7 and its appendix in the
  current manuscript for that proof.
- `provenance/canonical-manifest.json` and `canonical-README.txt`: unchanged
  records from the inherited source package. That historical manifest also
  lists an older PDF and older verification outputs that are intentionally not
  redistributed here; it is not this release's file inventory.
- `provenance/verified-source-sha256.json`: the 33 source hashes from the fresh
  replay. `check_inventory.py` verifies all of them against this checkout.
- `release-inventory.json`: sizes and SHA-256 hashes for every distributed file
  except the inventory itself. `check_inventory.py` checks this inventory,
  source provenance, and exported-report hashes without requiring python-flint.

No license has been selected for this repository.
