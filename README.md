# Correlated Forest Rounding for Unweighted Nash Social Welfare: verification code

This package accompanies the current manuscript. It contains the numerical
verification programs, exact inputs, certificates, provenance, and saved replay
reports used by the subadditive and tight-submodular analyses. It is not an
end-to-end allocation solver or a formal verification of all mathematical proofs.

## Contents and reproduction

- `tight-submodular/`: the complete, unchanged October 30, 2026 release for the
  tight-submodular analysis in Section 7 and Appendix F. Follow its `README.md`
  for dependency versions and full replay commands. Start in that directory;
  its `reproduce.py` calls `verify.py` and the additional constant audit.
- `verification/`: the unchanged exact subadditive concentration/calibration
  verifier and certificate used in Appendix E. Run
  `cd verification && python3 verify-subadditive.py` with Python 3.10 or later.
  No third-party dependency is needed for this verifier.
- `SHA256SUMS`: hashes of all files in this package except this checksum file.

The source references `verification/verify-subadditive.py` using the path retained
here. Other verifier names in Appendix F are relative to `tight-submodular/`.
The `paper/` directory inside that release is a historical mathematical-source
snapshot documenting the verifier, not the current manuscript. It remains
unchanged to preserve the release inventory and provenance.

## Verification status and scope

On November 2, 2026, all 75 files of the tight-submodular release were compared
byte-for-byte with the previously audited release, its inventory checker passed,
and all three subadditive files were compared with manuscript commit `65ebaa5`.
This packaging check did not rerun the numerical certificates. The saved reports
and the scope of their previous computations are documented in each directory.
A successful numerical replay does not establish the analytic reductions in the
paper; those remain separate proof obligations.

The older `supplement/` package concerned the superseded 1.9984 calibration. It
is not a dependency of this manuscript or of this code release and is retained
in the separate historical archives.

## Publishing

Upload the contents of this directory to the future GitHub repository. Replace
the existing `\ruilong` GitHub-link placeholder in `appendix.tex` with the
repository URL. No repository has been created or code uploaded by this task.
The paper source is distributed separately and contains no verification code.
