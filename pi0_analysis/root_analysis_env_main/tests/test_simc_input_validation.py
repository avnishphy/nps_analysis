"""Run after Hall C setup; arguments: PRODUCER FIXTURE_GENERATOR REPO_ROOT.

The caller compiles tests/simc_input_fixture.C and the producer. All artifacts
go to a new scratch directory. This exercises input safety, not physics fits.
"""
import os
from pathlib import Path
import subprocess
import sys
import tempfile


def main():
    producer, generator, repo = map(lambda p: str(Path(p).resolve()), sys.argv[1:])
    work = Path(tempfile.mkdtemp(prefix="simc_input_regression_"))
    base = os.environ.copy()
    base.update({
        "NPS_KIN": "KinC_x60_4b", "NPS_ACCEPTANCE_CUTS_CONFIG": str(Path(repo)/"config/acceptance_cuts.conf"),
        "NPS_EBEAM": "10.538", "NPS_Z_NPS_CM": "407", "NPS_NPS_THETA_DEG": "12.2",
        "NPS_SIMC_PRODUCTION_MODE": "nominal-only", "NPS_SMEARING_MODE": "off",
        "NPS_PIPELINE_RUN_ID": "input-regression", "NPS_PRODUCER_INPUT_SET_IDENTITY": "fixture",
    })
    for mode in ("valid", "wrong-type", "short-branch"):
        directory = work/mode
        subprocess.run([generator, str(directory), mode], check=True)
        env = base.copy()
        for kind, channel in (("EXCLUSIVE", "exclusive"), ("SIDIS", "sidis"), ("DELTA", "delta")):
            env[f"NPS_SIMC_{kind}_INPUT"] = str(directory/f"{channel}.root")
            env[f"NPS_SIMC_{kind}_NORMFAC"] = "1"
            env[f"NPS_SIMC_{kind}_NGEN"] = "40"
        output = directory/"nominal.root"
        env["NPS_SIMC_OUTPUT_FILE"] = str(output)
        result = subprocess.run([producer], env=env, cwd=repo, capture_output=True, text=True)
        (directory/"producer.log").write_text(result.stdout + result.stderr)
        expected = 0 if mode == "valid" else 7
        if result.returncode != expected:
            raise AssertionError(f"{mode}: status {result.returncode}, expected {expected}; see {directory}")
        if mode == "valid":
            subprocess.run([producer, "--validate-output", str(output), "nominal",
                            "input-regression", "fixture", "ignored", "120"], check=True)
            subprocess.run([generator, str(output), "check-output"], check=True)
        else:
            reason = "Cannot bind input branch 'Weight'" if mode == "wrong-type" else "incomplete input branch 'Weight'"
            if reason not in result.stderr:
                raise AssertionError(f"Missing branch diagnostic: {directory}")
        print(f"PASS {mode}: producer status={result.returncode}")
    print(f"Artifacts: {work}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main()
