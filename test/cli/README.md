# CLI tests

All CLI tests live under `test/cli/`. The implementation is only in `src/gwaslab_cli/` (the legacy `gwaslab.CLI` shim has been removed).

- **Smoke tests** (`test_cli_entry.py`): call `gwaslab_cli.main` in-process (version, `--help`, QC + output, `formatbook list`).
- **Subcommands & errors** (`test_cli_subcommands_errors.py`): `config` / `path` / `formatbook`, launcher edge cases, parser errors, plot/extract/liftover branches (including mocks where full runs are heavy).
- **Unified pipeline** (`test_cli_unified.py`, unittest): broad flag combinations on `test/raw` fixtures (QC, harmonize, formats, workflows).
- **Bash examples** (`test_cli_bash_examples.py`): run `examples/10_cli/*.sh` in an isolated temp tree with `../../test/...` paths preserved. The wrapper on `PATH` runs `python -m gwaslab_cli.main`, and `PYTHONPATH` includes this repo’s `src/`.

Run (from repo root):

```bash
PYTHONPATH=src pytest test/cli/ -m "not network"
```

Network-dependent examples (`06_extract.sh`) are marked `@pytest.mark.network`.
`08_utility.sh` is skipped in automated runs (GWAS Catalog download is slow); run it manually when needed.

```bash
PYTHONPATH=src pytest test/cli/ -m network
```

`04_format_conversion.sh` is skipped in CI-style runs because the bundled toy sumstats can produce VCF lines that fail `tabix` indexing; run that script manually when validating VCF + tabix.
