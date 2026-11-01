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
- [`VERIFICATION.md`](VERIFICATION.md): the current manuscript-to-code map,
  the numerical checks, and the separate analytic proof obligations.
- `SHA256SUMS`: hashes of all distributed files except this checksum file.

The source references `verification/verify-subadditive.py` using the path retained
here. The tables in Appendix F use names relative to
`tight-submodular/certificate/`; its replay entry points are in
`tight-submodular/`.
The `paper/` directory inside that release is a historical mathematical-source
snapshot documenting the verifier, not the current manuscript. It remains
unchanged to preserve the release inventory and provenance.

## Verification status and scope

On November 2, 2026, all 75 files of the tight-submodular release were compared
byte-for-byte with the previously audited release, its inventory checker passed,
and all three subadditive files were compared with manuscript commit `65ebaa5`.
Before this GitHub upload, a fresh full replay also passed: 4,912,892 bounded
boxes, all auxiliary gates, and the additional directed-constant audit. The
subadditive verifier passed all 282 intervals and reproduced the shipped
certificate byte-for-byte. See [REPLAY-2026-11-02.json](REPLAY-2026-11-02.json)
for this run summary and source hashes. The complete October 30 reports and
logs remain in the unchanged tight-submodular release.
A successful numerical replay does not establish the analytic reductions in the
paper; those remain separate proof obligations.

The older `supplement/` package concerned the superseded 1.9984 calibration. It
is not a dependency of this manuscript or of this code release and is retained
in the separate historical archives.

## Repository and manuscript

The code repository is
[Correlated_Forest_Rounding_Unweighted_NSW](https://github.com/ruilong-zhang/Correlated_Forest_Rounding_Unweighted_NSW).
The manuscript cites this repository in Appendix F and links the subadditive
verifier in Appendix E. The paper's LaTeX source is distributed separately;
`tight-submodular/paper/` is only the frozen historical snapshot described above.
No separate current manuscript or unverified later improvement is included here.
