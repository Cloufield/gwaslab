# Sumstats

The main GWASLab object for loading, QC, harmonization, filtering, plotting, and downstream analysis of GWAS summary statistics.

```python
import gwaslab as gl

mysumstats = gl.Sumstats("sumstats.txt.gz", fmt="plink2")
mysumstats.basic_check()
mysumstats.harmonize(ref_seq="ref.fa")
mysumstats.plot_mqq()
```

## API sections

| Section | Description |
|---------|-------------|
| [Core](core.md) | Construction, pipelines (`basic_check`, `harmonize`), harmonization helpers, I/O |
| [Fix](fix.md) | Standardization and QC (`fix_id`, `fix_chr`, `check_sanity`, …) |
| [Filter](filter.md) | Subsetting variants (`filter_value`, `filter_region`, `exclude_hla`, …) |
| [Plot](plot.md) | Visualization methods (`plot_mqq`, `plot_region`, …) |
| [Downstream](downstream.md) | Lead/novel extraction, LDSC, clumping, finemapping, PRS |

For tutorials see [Sumstats Object](../SumstatsObject.md), [QC and Filtering](../QC&Filtering.md), and [Visualization](../Visualization.md).
