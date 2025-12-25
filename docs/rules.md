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

## Admonitions

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
