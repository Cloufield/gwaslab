# Dtype transformation memory audit

Audit date: 2026-07-06. Follow-up to gwaslab issue [#222](https://github.com/Cloufield/gwaslab/issues/222) (categorical `.astype(str)` OOM on allele columns).

## Summary

| Tier | Sites found | Fixed | Accepted risk |
|------|-------------|-------|---------------|
| 1 — Critical | 2 pattern classes | CI-guarded + `categorical_safe_str` | 1 (`_fast_check_status` fixed-width path) |
| 2 — High | 4 | 3 | 1 (`check_status` global FASTA pad) |
| 3 — Medium | 5 | 4 | 0 (`_fix_allele` uses category-aware helpers) |
| 4 — Low | 4 | 1 (`util_in_simulate`) | 3 documented |

Regression coverage: [`test/test_categorical_astype_memory.py`](../../test/test_categorical_astype_memory.py) — 19 tests (11 tracemalloc + 3 helper unit + 5 static-scan guards).

---

## Risk taxonomy

### Tier 1 — Critical

| Pattern | Memory formula | Trigger |
|---------|----------------|---------|
| `Categorical.astype(str)` | `n_rows × max_category_length` | One long SV allele in category index |
| `.to_numpy().astype(f'<U{max_len}')` on full column | `n_rows × max_len × 4` bytes | Same; used in harmonize ref-check |

### Tier 2 — High

| Pattern | Memory formula | Trigger |
|---------|----------------|---------|
| Full-column `categorical_safe_str` when subset suffices | `2 × n_rows × avg_len` temporary | Flip RC on few rows but materialize all EA/NEA |
| `np.pad(record, max_padding_needed)` | `genome_len + max_allele_len` | Global max allele length in batch |

### Tier 3 — Medium

| Pattern | Notes |
|---------|-------|
| `.str.len()` / `.str.upper()` on categorical without safe convert | Usually OK via categorical accessor; fragile across pandas versions |
| `pd.Categorical(materialized_series, ...)` after full string column | Temporary 2× peak |

---

## Per-site findings

### Mitigated (baseline)

| File | Function | Fix |
|------|----------|-----|
| `qc/qc_check_datatype.py` | `categorical_safe_str` | Project default; `astype("string")` |
| `hm/hm_casting.py` | `_merge_mold_with_sumstats_by_chrpos` | Uses `categorical_safe_str` for _ASET |
| `util/util_ex_ldsc.py` | `_munge_sumstats` | Safe allele concat |
| `util/util_in_filter_value.py` | `_filter_values`, `_filter_palindromic` | Safe paths + tracemalloc |
| `extension/.../validate_ssf.py` | `_validate_dataframe` | Safe allele regex checks |
| `hm/hm_harmonize_sumstats.py` | `is_palindromic` | Safe before `to_numpy()` |
| `qc/qc_fix_sumstats.py` | Most ID/allele fix paths | `astype("string")` / `categorical_safe_str` |

### Fixed in this audit

| File | Function | Change | Tier |
|------|----------|--------|------|
| `io/io_to_formats.py` | `_check_indel` | Wrap EA/NEA with `categorical_safe_str` before `.str.len()` | 3 |
| `io/io_to_formats.py` | `_adjust_position` | Use safe strings for len/slice ops | 3 |
| `hm/hm_assign_rsid.py` | `_expand_multiallelic_fast` | Safe convert before `.str.contains` / `.str.split` | 3 |
| `util/util_in_simulate.py` | VCF extract | `np.asarray(..., dtype=object)` instead of `.astype(str)` | 4 |
| `qc/qc_fix_sumstats.py` | `_flip_allele_stats` RC branch | Use `.cat.categories` + `add_categories`; no full-column materialization | 2 |
| `qc/qc_fix_sumstats.py` | `_fix_allele`, `normalizae_chunk` | `categorical_str_len` / `categorical_str_upper` / `categorical_str_contains` | 3 |
| `util/util_in_filter_value.py` | `_filter_indel`, `_filter_snp` | `categorical_str_len` | 3 |
| `hm/hm_infer_with_af.py` | strand inference | `categorical_str_len` | 3 |
| `hm/hm_harmonize_sumstats.py` | `check_status`, `_fast_check_status` | `categorical_str_len` for partition/max len | 3 |
| `qc/qc_check_datatype.py` | `categorical_str_len`, `categorical_str_upper`, `categorical_str_contains` | Category-index helpers | 3 |

### Accepted risk (documented)

| File | Function | Rationale |
|------|----------|-----------|
| `hm/hm_harmonize_sumstats.py` | `_fast_check_status` | Intentional fixed-width numpy path for vectorized ref alignment. **Mitigation:** `check_status` partitions at `max_len=4` so SNPs stay in small-padded batch; long alleles isolated to second pass. Allowlisted in CI (`TO_NUMPY_UNICODE_ALLOWLIST`). |
| `hm/hm_harmonize_sumstats.py` | `check_status` L850 | Pads concatenated FASTA buffer by global max allele length in batch. Cost is O(genome), not O(n_rows × max_len) for the buffer itself. |
| `hm/hm_casting.py` | `_ASET` construction | Full-column `categorical_safe_str` required for merge keys. Allowlisted (`FULL_COLUMN_CATEGORICAL_SAFE_STR_ALLOWLIST`). |

### Low priority (no change)

| File | Pattern | Notes |
|------|---------|-------|
| `viz/viz_plot_track.py` | `.index.astype(str)` | Small index |
| `info/g_Sumstats_summary.py` | `columns.astype(str)` | Metadata only |
| `io/io_fasta.py` | `np.array([r], dtype=f'<U{r_len}')` | Per-chromosome, not per-variant |
| `extension/ldsc/ldsc_sumstats.py` | `alleles.apply(lambda ...)` | Small fixed allele-key set |

---

## CI guards (`TestDtypeMemoryPatterns`)

| Test | Pattern | Action on violation |
|------|---------|---------------------|
| `test_no_unsafe_astype_str_on_allele_columns` | `[EA\|NEA\|REF\|ALT].astype(str)` | FAIL |
| `test_prefer_pandas_string_in_sumstats_modules` | `.astype(str)` in sumstats modules | FAIL (with allowlist) |
| `test_to_numpy_unicode_astype_is_allowlisted` | `.to_numpy().astype('<U...')` in hm/qc/io | FAIL if not in allowlist |
| `test_allele_str_ops_use_safe_conversion` | `[EA\|NEA\|REF\|ALT].str.*` without nearby safe helper | FAIL (allowlist empty after migration) |
| `test_no_full_column_categorical_safe_str_on_alleles` | `categorical_safe_str(sumstats['EA'])` in qc/hm | FAIL (with allowlist) |

---

## tracemalloc regression matrix

Fixture: 100,000 rows, categorical EA/NEA, one 662 bp allele in category index. Threshold: 100 MB peak.

| Function | File | Status |
|----------|------|--------|
| `_merge_mold_with_sumstats_by_chrpos` | hm_casting | PASS |
| `_munge_sumstats` | util_ex_ldsc | PASS |
| `_filter_values` | util_in_filter_value | PASS |
| `_validate_dataframe` | validate_ssf | PASS |
| `is_palindromic` | hm_harmonize_sumstats | PASS |
| `_filter_palindromic` | util_in_filter_value | PASS |
| `check_status` | hm_harmonize_sumstats | PASS |
| `_fix_allele` | qc_fix_sumstats | PASS |
| `_flip_allele_stats` | qc_fix_sumstats | PASS |
| `_check_indel` | io_to_formats | PASS |
| `_expand_multiallelic_fast` | hm_assign_rsid | PASS |
| `_filter_indel` | util_in_filter_value | PASS |
| `_filter_snp` | util_in_filter_value | PASS |

Repro notebook: [`examples/bug/replicate_issue_222_astype_str_memory.ipynb`](../../examples/bug/replicate_issue_222_astype_str_memory.ipynb).

---

## Recommendations for future code

1. **Store** alleles as `category` when possible.
2. **Convert** with `categorical_safe_str()` only at string-op boundaries (concat, `.str`, export).
3. **Never** use `.astype(str)` on categorical EA/NEA/REF/ALT.
4. **Prefer** `.cat.categories` / `add_categories` or `categorical_str_len` / `categorical_str_upper` over full-column string materialization.
5. **Add tracemalloc test** when introducing new allele string conversion in hot paths.
