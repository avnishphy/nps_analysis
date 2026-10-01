# Provenance manifests

Tracked manifests describe immutable inputs and the relationship between
results and source/configuration. Generated campaigns must record the baseline
commit, implementation commit, configuration hashes, input identities,
software versions, commands, RNG streams, and output checksums.

`baseline_copy_manifest.tsv` maps every migrated file to its baseline Git blob.
`baseline_copy_manifest.SHA256SUMS` protects that mapping. The full ignored and
untracked `MAIN` freeze manifests are preserved under
`recovery/main_checkpoint_20261001/` and intentionally remain outside Git.
