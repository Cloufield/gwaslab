# Documentation Style Guide

This document outlines the style and formatting rules for GWASLab documentation.

!!! note "Documentation System Migration"
    The documentation system is transitioning from mkdocs-material to zensical. From now on, use zensical for all new documentation work.

## Code Examples

- **Code block types**: Use Python code blocks for GWASLab Python API code; use bash code blocks for GWASLab CLI commands.
- **Variable naming**: Always use `mysumstats` as the standard variable name for Sumstats objects in sample code. Use `mysumstats.function()` instead of `.function()` or other variable names.

## Formatting

- **Lists**: There should be an empty line before any lists
- **Reserved headers**: Column names and reserved headers (like `SNPID`, `CHR`, `POS`, `EA`, `NEA`, `BETA`, `SE`, `P`, etc.) should be in **bold** format: `**SNPID**`, `**CHR**`, etc.
- **Table headers**: Use `DataType` (camelCase, no space) consistently in all parameter/option tables, not `Type` or `Data Type`

## Docstrings

All Python docstrings must use **NumPy (numpydoc)** format. This matches [mkdocstrings](https://mkdocstrings.github.io/) (`docstring_style: numpy`) and drives the [API Reference](api/index.md).

- **Section headers**: `Parameters`, `Returns`, `Notes`, `Raises`, `Examples` — no trailing colon
- **Underlines**: match the header length (`Parameters` → ten dashes, `Returns` → seven, etc.)
- **Parameters with defaults**: `name : bool, default False` (space after `default`, not `default=`)
- **Parameters without defaults**: `name : str, optional`
- **Descriptions**: indented four spaces under the parameter line

Canonical template:

```python
"""
One-line summary.

Extended description if needed.

Parameters
----------
name : bool, default False
    What it does.
other : str, optional
    No default; describe behavior.

Returns
-------
Sumstats or None
    When ``inplace=True``, returns None.

Notes
-----
Additional context.

Examples
--------
>>> import gwaslab as gl
>>> mysumstats = gl.Sumstats("sumstats.txt.gz", fmt="plink2")
>>> mysumstats.basic_check()
```

Plot function docstrings are generated from the visualization parameter registry via `gwaslab.info.g_numpy_doc` and `gwaslab.viz.viz_aux_doc`. Sumstats method wrappers use `@add_doc` in `gwaslab.info.g_object_helper` (adapted from implementation docstrings). For bulk docstring cleanup, use `gwaslab.info.g_numpy_doc.normalize_sections`.

See also `src/gwaslab/dev_principles.txt` (Documentation section).


!!! note "Note on Admonitions"
    The documentation system is transitioning from mkdocs-material to zensical. Use appropriate admonition blocks for different types of content:

- **Examples**: Use `!!! example` blocks for code examples and usage demonstrations
- **Citations/References**: Use `!!! quote` blocks for citations, references, and quoted material
- **General explanations**: Use `!!! note` blocks for general explanations, tips, and additional information
- **Info**: Use `!!! info` blocks for informational content and announcements
- **Warnings**: Use `!!! warning` blocks for important warnings and cautions

## Examples

### Code Block Format
````python
# Correct
```python
mysumstats.plot_mqq()
```

# Incorrect  
```
mysumstats.plot_mqq()
```
````

### Variable Naming
```
# Correct
mysumstats.basic_check()
mysumstats.plot_mqq()

# Incorrect
sumstats.basic_check()
.plot_mqq()
.function()
```python

### Reserved Headers
```
# Correct
The **SNPID** column contains variant identifiers.
Use **CHR** and **POS** for genomic coordinates.

# Incorrect
The SNPID column contains variant identifiers.
Use CHR and POS for genomic coordinates.
```python

### Table Format
```
# Correct
| Parameter | DataType | Description | Default |

# Incorrect
| Parameter | Type | Description | Default |
| Parameter | Data Type | Description | Default |
```python
