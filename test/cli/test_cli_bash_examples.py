"""Run `examples/10_cli/*.sh` the same way users do (via `bash` + `gwaslab` on PATH)."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

from gwaslab_cli.main import main

REPO_ROOT = Path(__file__).resolve().parents[2]
EXAMPLES_CLI_DIR = REPO_ROOT / "examples" / "10_cli"

# Mirrors examples shipped with the repo (no network).
BASH_EXAMPLES_OFFLINE = [
    "01_basic_qc.sh",
    "02_harmonize.sh",
    "03_liftover.sh",
    "04_format_conversion.sh",
    "05_plot.sh",
    "07_pipeline.sh",
]

# Needs GWAS Catalog / external services.
BASH_EXAMPLES_NETWORK = [
    "06_extract.sh",
    "08_utility.sh",
]


def test_cli_example_scripts_present() -> None:
    """Ensure every documented example script exists under examples/10_cli/."""
    for script in BASH_EXAMPLES_OFFLINE + BASH_EXAMPLES_NETWORK:
        path = EXAMPLES_CLI_DIR / script
        assert path.is_file(), f"Missing example script: {path}"


@pytest.mark.parametrize("script", BASH_EXAMPLES_OFFLINE)
def test_bash_example_offline(
    script: str,
    example_cli_workdir: Path,
    gwaslab_env: dict[str, str],
) -> None:
    """Offline examples should complete with exit code 0."""
    if script == "04_format_conversion.sh":
        pytest.skip(
            "04_format_conversion.sh: VCF + --tabix step can fail on toy sumstats (unsorted "
            "coordinates / tabix indexing). Other steps are covered by formatbook + conversion tests."
        )

    proc = subprocess.run(
        ["bash", str(example_cli_workdir / script)],
        cwd=str(example_cli_workdir),
        env=gwaslab_env,
        capture_output=True,
        text=True,
        timeout=900,
    )
    assert proc.returncode == 0, (
        f"{script} failed\n--- stdout ---\n{proc.stdout}\n--- stderr ---\n{proc.stderr}"
    )


@pytest.mark.network
@pytest.mark.parametrize("script", BASH_EXAMPLES_NETWORK)
def test_bash_example_network(
    script: str,
    example_cli_workdir: Path,
    gwaslab_env: dict[str, str],
) -> None:
    """Examples that call remote APIs or `formatbook update` / downloads."""
    if script == "08_utility.sh":
        pytest.skip(
            "08_utility.sh: GWAS Catalog download (GCST90270926) is slow and flaky in CI; "
            "run manually with: pytest test/cli/test_cli_bash_examples.py -m network -k 08_utility"
        )

    proc = subprocess.run(
        ["bash", str(example_cli_workdir / script)],
        cwd=str(example_cli_workdir),
        env=gwaslab_env,
        capture_output=True,
        text=True,
        timeout=900,
    )
    assert proc.returncode == 0, (
        f"{script} failed\n--- stdout ---\n{proc.stdout}\n--- stderr ---\n{proc.stderr}"
    )


def test_gwaslab_module_invocation_matches_entrypoint(gwaslab_env: dict[str, str]) -> None:
    """`python -m gwaslab_cli.main version` should work (used by bash wrapper)."""
    proc = subprocess.run(
        [sys.executable, "-m", "gwaslab_cli.main", "version"],
        capture_output=True,
        text=True,
        timeout=60,
        env=gwaslab_env,
    )
    assert proc.returncode == 0, proc.stderr


def test_main_version_matches_subprocess(
    monkeypatch: pytest.MonkeyPatch,
    gwaslab_env: dict[str, str],
) -> None:
    """In-process `main(["version"])` should match subprocess module call."""
    monkeypatch.setattr(sys, "argv", ["gwaslab", "version"])
    main(["version"])
    p = subprocess.run(
        [sys.executable, "-m", "gwaslab_cli.main", "version"],
        capture_output=True,
        text=True,
        timeout=60,
        env=gwaslab_env,
    )
    assert p.returncode == 0
