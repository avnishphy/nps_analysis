# External and excluded input register

Status: initial migration inventory; active campaign selection remains
unresolved.

| Input | Role | Identity/provenance | Policy |
|---|---|---|---|
| `docs/pi0_uncertainty_audit_2026-10-01.md` in the older `root_analysis_env` workspace | Scientific audit | External path and SHA-256 are recorded in `read_only_reference_SHA256SUMS` | Read-only; never reconstructed or edited here |
| `/group/nps/singhav/setup.csh` plus module stack | Software environment | ROOT 6.30.04 and captured versions in `software_environment.txt` | Source before ROOT/build commands |
| Frozen `MAIN/output/` tree | Historical generated artifacts | Full file hashes are in ignored checkpoint recovery manifests | Read-only; do not use as an implicit production destination; provenance must be established before any baseline input is selected |
| Hall C replay/production files named by analysis launchers | Event input | Exact files, segments, sizes, hashes/catalog IDs, and coverage are not yet inventoried | Read-only; campaign manifest required before execution |
| SIMC/Geant4 products under external `/volatile` or SIMC work areas | Simulation input | Existing source contains several absolute candidates; active versions are not yet selected | Read-only; require channel/run/config manifest and generator/stopping provenance |
| Primary literature/reference PDFs retained only in the frozen efficiency research tree | Method evidence | Covered by full `MAIN` freeze hashes; not active numerical inputs | Keep read-only; cite exact source/section when used |

Small active configuration/calibration inputs migrated into `config/` are
checksum-protected by `config_reference_SHA256SUMS`. Generated outputs and raw
efficiency caches were not copied into the committed baseline snapshot; their
frozen identities remain recoverable through the checkpoint manifests.
