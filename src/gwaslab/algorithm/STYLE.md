# Algorithm module style guide

Apply to every function under `gwaslab/algorithm/`.

## Function style

- `from __future__ import annotations` in every module
- Full type hints; inputs coerced with `np.asarray(..., dtype=float)` at boundaries
- No `Sumstats`, `Log`, matplotlib, file I/O, or `@with_logging`
- Public names: `snake_case`; private helpers: leading `_`
- Constants: `UPPER_SNAKE` at module level when shared

## Docstrings (NumPy format)

1. One-line summary (imperative)
2. Extended summary (optional)
3. Parameters
4. Returns
5. Notes (orchestration pointers, GWASLab-only heuristics)
6. References (required for published methods)

Example: see `core/genomic_control.py`.

## Citations

- `References` section with Author (Year). Title. Journal, volume(issue), pages.
- GWASLab-only algorithms: document in `Notes`; omit empty References

## Comments

- Module docstring: purpose + orchestration layer pointer
- Inline comments only for non-obvious math
- No date stamps, TODO/FIXME, or `# ----` banners

## PR checklist

- [ ] Type hints on all public functions
- [ ] NumPy docstring with References when applicable
- [ ] No forbidden imports (`g_Log`, `g_Sumstats`, `matplotlib`)
- [ ] Unit test in `test/algorithm/`
- [ ] Existing integration tests still pass
