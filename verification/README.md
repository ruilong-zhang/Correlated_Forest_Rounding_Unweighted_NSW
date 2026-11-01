# Subadditive concentration and calibration certificate

Run `python3 verify-subadditive.py` from this directory (Python 3.10 or later; no third-party packages). It writes `certificate.json` and reports whether all exact comparisons pass.

The source is the full note's Appendix D verifier, with only its introductory description updated. Every decisive comparison uses rational arithmetic. The manuscript proves the reductions that make the 282 interval comparisons cover all real means; these checks are not samples of the mean parameter. The script also encloses the final approximation constant strictly between 35.59743 and 35.59744.

This certificate verifies the scalar concentration/calibration inequalities. It is not a formal verification of the configuration optimization, contention, or forest-rounding algorithms.
