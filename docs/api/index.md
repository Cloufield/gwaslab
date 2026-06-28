# API Reference

This section is auto-generated from **NumPy (numpydoc)** docstrings in the GWASLab source code using [mkdocstrings](https://mkdocstrings.github.io/). Parameter tables use a consistent layout (`docstring_section_style: table`). See [Documentation Style Guide — Docstrings](../rules.md#docstrings) for the canonical format.

## Quick start

```python
import gwaslab as gl

mysumstats = gl.Sumstats("sumstats.txt.gz", fmt="plink2")
mysumstats.basic_check()
mysumstats.plot_mqq()
```

## Pages

| Page | Contents |
|------|----------|
| [Package functions](package.md) | Top-level `gl.*` re-exports (non-plot) |
| [Reference & utilities](reference.md) | `gl.download_ref`, `gl.get_path`, `gl.get_power`, … |
| [I/O helpers](io.md) | `gl.read_ldsc`, `gl.load_pickle`, tabular readers |
| [Plotting](plotting.md) | Top-level `gl.plot_*`, `gl.compare_effect`, `gl.scatter` — user guides: [Lead Overlap](../LeadOverlapPlot.md), [Sankey](../SankeyPlot.md) |
| [Sumstats → Core](sumstats/core.md) | `gl.Sumstats` loading, pipelines, harmonization, I/O |
| [Sumstats → Fix](sumstats/fix.md) | `fix_id`, `fix_chr`, `check_sanity`, … |
| [Sumstats → Filter](sumstats/filter.md) | `filter_value`, `filter_region`, `exclude_hla`, … |
| [Sumstats → Plot](sumstats/plot.md) | `mysumstats.plot_mqq`, `plot_region`, … |
| [Sumstats → Downstream](sumstats/downstream.md) | `get_lead`, LDSC, clumping, finemapping, PRS |

For tutorials and workflow examples, see [Tutorial](../tutorial_v4.md), [QC and Filtering](../QC&Filtering.md), [Harmonization](../Harmonization.md), and [Visualization](../Visualization.md).

## Build locally

```bash
pip install -e ".[docs]"
pip install zensical
zensical serve   # live preview
# or:
zensical build --clean
```

!!! note "Generated docs"
    Parameter lists for plot functions reflect the visualization parameter registry and match `help(gl.plot_miami2)` / `help(mysumstats.plot_region)`.
