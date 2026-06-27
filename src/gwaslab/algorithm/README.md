# GWASLab `algorithm/` package

Pure calculation layer for GWAS summary-statistics math. No `Sumstats`, no logging, no file I/O.

## Layout

| Subpackage | Role |
|------------|------|
| `core/` | Stat conversions, lambda GC, power, HWE, h2 |
| `heterogeneity/` | Meta-analysis, Cochran Q, jackknife correlation |
| `finemapping/` | Wakefield ABF, winner's curse |
| `allele/` | Reverse-complement, palindromic tests |
| `leads/` | Distance-window lead clustering |
| `density/` | Brisbane / signal density |
| `population/` | F_ST and related population metrics |

## Boundaries

- **Orchestration:** `util_in_*`, `qc/`, `hm/`, `viz/` — DataFrames, logging, Sumstats methods
- **Large tools:** `extension/` — LDSC, MR, MTAG, PRS-CS, MultiSuSiE (unchanged)

## Style

All public functions follow [`STYLE.md`](STYLE.md): NumPy docstrings, type hints, literature `References` where applicable.
