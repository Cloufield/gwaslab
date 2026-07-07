# Simulate GWAS summary statistics

!!! info "In development"

    `simulate_sumstats_region` and `simulate_sumstats_global` are under active development. These pages are not yet listed in the site navigation. Runnable workflow: [Simulate sumstats workflow](utility_simulate_sumstats.md). Development notebook: `examples/_development/simulation/utility_simulate_sumstats.ipynb`.

GWASLab can generate **simulated GWAS summary statistics** (BETA, SE, Z, P) from a **reference-panel VCF** using an efficient LD model based on standardized genotypes (`X^T X`), without building full LD matrices in memory.

**Key features:**

- Region-scoped or genome-wide simulation
- Quantitative or case–control traits
- Sparse or polygenic genetic architecture
- Optional MAF-dependent effect sizes (`alpha`)
- Realism knobs: sample-size dropout, imputation INFO, λ<sub>GC</sub>, population stratification

## Functions

| Function | Scope | Use when |
|----------|-------|----------|
| [`simulate_sumstats_region`](#simulate_sumstats_region) | One genomic region | Quick demos, regional plots, method tests |
| [`simulate_sumstats_global`](#simulate_sumstats_global) | One or more chromosomes | Genome-wide workflows; applies global heritability (`h2`) calibration |

Both return `(Sumstats, causal_snp_ids)`:

```python
from gwaslab import simulate_sumstats_region, simulate_sumstats_global

sumstats, causals = simulate_sumstats_region(
    vcf_path="path/to/panel.vcf.gz",
    region=("7", 126753550, 127753550),
    n=10_000,
    n_causal=3,
    seed=42,
)
```

The `Sumstats` object behaves like any other GWASLab sumstats object (`basic_check`, `plot_mqq`, export, etc.).

## Prerequisites

| Requirement | Notes |
|-------------|-------|
| Reference VCF/BCF | Genotypes in `FORMAT/GT`; allele frequency in INFO (e.g. `AF`) |
| Tabix index | `.tbi` beside the VCF for region queries |
| Region tuple | `(chrom, start, end)` with 1-based coordinates for region mode |

Bundled demo panel (also used in tests and [tutorial v4](tutorial_v4.md)):

```
test/ref/1kg_eas_hg19.chr7_126253550_128253550.vcf.gz
```

Or from sample data / registry:

```python
vcf_path = "gwaslab-sample-data/1kg_eas_hg19.chr7_126253550_128253550.vcf.gz"
# vcf_path = gl.get_path("1kg_eas_hg19")  # full genome; pass region= or chromosomes=
```

## Model (brief)

Simulated Z-scores follow:

```
z = μ_causal + b_strat + ε
```

- **μ_causal** — causal signal propagated through LD (`sqrt(N_eff) * (1/n) X^T (X β)`)
- **ε** — noise scaled by cryptic relatedness (`lambda_gc`)
- **b_strat** — optional stratification bias (`sigma_strat`)

Effect sizes are drawn on a standardized scale; summary columns are derived from Z and effective sample size. See docstrings in `gwaslab.util.util_in_simulate` for the full mathematical description.

### Model assumptions

| Concept | Source | Notes |
|---------|--------|-------|
| LD structure | Reference-panel genotypes (`n_ref` samples in the VCF) | Standardized `X`; noise uses `(1/√n_ref) X^T u` |
| Association strength | User GWAS `N` or `N_eff` | Causal mean scales with `√N_eff`; SE = `1/√N_eff` |
| Heritability (`h2`) | Global mode (required); region mode (optional `h2=`) | Rescales `BETA_TRUE` so `(1/n_ref) ‖Xβ‖²` matches target |
| Small panels | Same formulas | Fewer reference samples → noisier LD; prefer larger panels when possible |

When `verbose=True`, simulations log `n_ref`, variant counts, `n_matched`, causal count, and (global) `V_g_raw`, `scale_factor`, and achieved `h2`.

**Binary traits:** default `trait_model="linear"` uses case–control `N_eff` with linear Z. `trait_model="liability"` applies a simplified probit threshold scale to `μ_causal` (documented limitation; not a full logistic GWAS).

**Global blocks:** LD noise is simulated independently within each `window_bp` block; block boundaries are not smoothed unless you use a smaller window or future overlap options.

## Output columns

| Column | Description |
|--------|-------------|
| `CHR`, `POS`, `EA`, `NEA`, `SNPID`, `EAF` | Variant identity and allele frequency |
| `Z`, `BETA`, `SE`, `P`, `MLOG10P` | Association statistics |
| `N`, `N_EFF` | Reported and effective sample size |
| `INFO` | Imputation quality (simulated or from VCF) |
| `IS_CAUSAL`, `BETA_TRUE` | Ground-truth labels for evaluation |
| `N_CASE`, `N_CONTROL` | Present for `trait="binary"` |

## `simulate_sumstats_region`

Simulate all variants in a region (optionally thinned).

### Trait and sample size

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `region` | `tuple` | required | `(chr, start, end)` |
| `trait` | `str` | `"quant"` | `"quant"` or `"binary"` |
| `n` | `int` | `300_000` | Sample size (quantitative) |
| `n_case` | `int` | `None` | Cases (binary) |
| `n_ctrl` | `int` | `None` | Controls (binary) |

### Genetic architecture

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `mode` | `str` | `"sparse"` | `"sparse"` (fixed `n_causal`, `0` allowed) or `"polygenic"` (Bernoulli `pi` per variant) |
| `n_causal` | `int` | `1` | Causal count when `mode="sparse"` |
| `pi` | `float` | `2e-3` | Causal probability when `mode="polygenic"` |
| `h2` | `float` or `None` | `None` | Optional local heritability rescale (region mode) |
| `trait_model` | `str` | `"linear"` | `"linear"` or `"liability"` (binary only) |
| `effect_sd` | `float` | `0.05` | Effect SD when `alpha=0` |
| `alpha` | `float` | `0` | MAF-dependent scaling; `>0` gives rarer variants larger effects |

### Filtering and realism

Neutral defaults (no artificial inflation). For demo-style plots, pass explicit realism kwargs (see [workflow](utility_simulate_sumstats.md#realism-knobs)).

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `maf_min`, `maf_max` | `float` | `0.01`, `0.5` | Variant MAF filters |
| `thin` | `float` or `None` | `None` | Fraction of variants to keep after filtering; `None` = all |
| `lambda_gc` | `float` | `1.0` | Cryptic relatedness / genomic-control inflation |
| `sigma_strat` | `float` | `0.0` | Stratification bias |
| `n_drop_rate` | `float` | `0.0` | Fraction of SNPs with reduced `N` |
| `use_info` | `bool` | `True` | Attenuate `N_EFF` by INFO |
| `seed` | `int` | `1` | Random seed |

### Metadata

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `study` | `str` | `"Simulated_Study"` | Study label |
| `trait_name` | `str` | `"Simulated_Trait"` | Trait label |
| `build` | `str` | `"38"` | Genome build metadata |
| `verbose` | `bool` | `True` | Log progress |

## `simulate_sumstats_global`

Genome-wide (or per-chromosome) simulation with block-wise processing and **global heritability calibration** (`h2`).

Additional options beyond region mode:

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `chromosomes` | `list` | all in VCF | e.g. `["7"]` or `[1, 2, 3]`; autodetected from VCF header when `None` |
| `h2` | `float` | `0.1` | Target heritability on standardized scale (all chromosomes) |
| `window_bp` | `int` | `1_000_000` | Block size for LD-panel windows |
| `null_mode` | `str` | auto | `"iid_far"` for sparse (fast), `"ld_panel"` for polygenic / strict LD |
| `ld_window_bp` | `int` | `window_bp` | Radius around each causal for LD panel when `null_mode="iid_far"` |
| `thin`, `lambda_gc`, `sigma_strat` | — | same as region | Applied in Phase 3 Z simulation |
| `trait_model` | `str` | `"linear"` | Binary trait scaling (see above) |
| `block_overlap_bp` | `int` | `0` | Reserved; block overlap not implemented yet |

!!! note "Sparse fast path (`null_mode`)"
    When `mode="sparse"`, the default `null_mode="iid_far"` assigns i.i.d. null Z-scores outside causal LD windows and uses the reference panel only near hits. This makes full-genome sparse demos practical (minutes instead of hours). Use `null_mode="ld_panel"` for LD-correlated null noise genome-wide.

!!! warning "Runtime"
    Full-genome `ld_panel` simulation on large reference panels can take substantial time and memory. Sparse simulations default to `iid_far`; start with a single chromosome for strict LD-panel tests.

## Visualization

The returned `Sumstats` object supports the usual downstream steps:

- `basic_check()` — validate columns before plotting or export
- `plot_mqq()` — Manhattan (`mode="m"`), QQ (`mode="qq"`), or regional (`mode="r"`) plots

Suggested workflow:

| Simulation | Plot | Key options |
|------------|------|-------------|
| `simulate_sumstats_region` | Regional | `region=`, `vcf_path=`, `pinpoint=causal_ids` |
| `simulate_sumstats_global` | Manhattan | `mode="m"`, `highlight=causals` |

Runnable examples with figures: [Simulate sumstats workflow](utility_simulate_sumstats.md#visualize-simulated-results).

## Examples

!!! example "Region, sparse architecture"
    ```python
    import gwaslab as gl

    vcf = "test/ref/1kg_eas_hg19.chr7_126253550_128253550.vcf.gz"
    sumstats, causals = gl.simulate_sumstats_region(
        vcf_path=vcf,
        region=("7", 126753550, 127753550),
        n=10_000,
        n_causal=3,
        mode="sparse",
        seed=42,
        build="19",
        verbose=False,
    )
    print(len(sumstats.data), causals)
    ```

!!! example "Global, full genome"
    ```python
    vcf_gw = gl.get_path("1kg_eas_hg19")
    sumstats, causals = gl.simulate_sumstats_global(
        vcf_path=vcf_gw,
        n=10_000,
        mode="sparse",
        n_causal=3,
        h2=0.01,
        seed=42,
        build="19",
        verbose=False,
    )
    ```

Runnable step-by-step examples: [Simulate sumstats workflow](utility_simulate_sumstats.md).

Development notebook: `examples/_development/simulation/utility_simulate_sumstats.ipynb`.
