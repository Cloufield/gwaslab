# GWASLab License Summary

## Project License

GWASLab is licensed under **GPL-3.0-only** (GNU General Public License v3.0).

See the root `LICENSE` file for the full GPL-3.0 license text.

### License History

**Note**: Prior to version 3.4.39, GWASLab was licensed under the MIT License. The license was changed to GPL-3.0-only starting with version 3.4.39. The previous MIT license text is preserved below for reference:

```python
MIT License

Copyright (c) 2022 Cloufield

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

**For versions < 3.4.39**: The MIT license applies.  
**For versions >= 3.4.39**: The GPL-3.0-only license applies.

## Third-Party Dependencies and Their Licenses

This document summarizes the licenses of all third-party packages and extensions used in GWASLab. All listed dependencies are compatible with GPL-3.0.

### Core Dependencies

| Package | License | Compatibility with GPL-3.0 | Notes |
|---------|---------|---------------------------|-------|
| **pandas** (>=1.3, !=1.5) | BSD-3-Clause | ✅ Compatible | Permissive license, can be used with GPL |
| **numpy** (>=1.21.2, <2) | BSD-3-Clause | ✅ Compatible | Permissive license, can be used with GPL |
| **matplotlib** (>=3.8, <3.9) | PSF-based (Python Software Foundation License) | ✅ Compatible | PSF license is GPL-compatible |
| **seaborn** (>=0.12) | BSD-3-Clause | ✅ Compatible | Permissive license, can be used with GPL |
| **scipy** (>=1.12) | BSD-3-Clause | ✅ Compatible | Permissive license, can be used with GPL |
| **pysam** (==0.22.1) | MIT License | ✅ Compatible | MIT is GPL-compatible |
| **adjustText** (>=0.7.3, <=0.8) | MIT License | ✅ Compatible | MIT is GPL-compatible |
| **scikit-allel** (>=1.3.5) | MIT License | ✅ Compatible | MIT is GPL-compatible |
| **h5py** (>=3.10.0) | BSD-3-Clause | ✅ Compatible | Permissive license, can be used with GPL |
| **pyarrow** | Apache-2.0 | ✅ Compatible | Apache-2.0 is GPL-3.0 compatible |
| **polars** (>=1.27.0) | MIT License | ✅ Compatible | MIT is GPL-compatible |
| **sumstats-liftover** (==1.1.0) | MIT License | ✅ Compatible | MIT is GPL-compatible |
| **requests** | Apache-2.0 | ✅ Compatible | Used for GWAS Catalog API access |

### Build System Dependencies

| Package | License | Compatibility with GPL-3.0 |
|---------|---------|---------------------------|
| **setuptools** (>=68) | MIT License | ✅ Compatible |
| **wheel** | MIT License | ✅ Compatible |

### Included Extension Code

GWASLab includes implementations of several analysis tools in the `extension/` directory:

| Component | Source | License | Notes |
|-----------|--------|---------|-------|
| **LDSC (LD Score Regression)** | Included in `gwaslab/extension/ldsc/` | GPL-3.0 | Based on LD Score regression methodology. Code included in GWASLab repository. |
| **PRS-CS (Polygenic Risk Score - Continuous Shrinkage)** | Included in `gwaslab/extension/prscs/` | MIT | PRS-CS implementation included in GWASLab. Original reference: Ge et al., Nature Communications, 10:1776, 2019. |
| **MultiSuSiE (Multi-ancestry SuSiE)** | Included in `gwaslab/extension/multisusie/` | GPL-3.0 | MultiSuSiE implementation included in GWASLab. Original reference: Rossen et al., Nature Genetics, 2026. Source: https://github.com/jordanero/MultiSuSiE |

### Optional External Tools (Runtime Dependencies)

The following external command-line tools are used by various `util_ex_*` functions but are **not Python packages** and are **not required** for basic GWASLab functionality. They are optional runtime dependencies:

| Tool | Used By | License | Notes |
|------|---------|---------|-------|
| **R** | All `util_ex_run_*` functions | GPL-2.0 or GPL-3.0 | R statistical computing environment |
| **susieR** | `util_ex_run_susie`, `util_ex_run_coloc`, `util_ex_run_mesusie` | BSD-style (see note) | Sum of Single Effects model for fine-mapping. Copyright (c) 2017-2022, Gao Wang, Peter Carbonetto, Yuxin Zou, Kaiqian Zhang, Matthew Stephens |
| **coloc** | `util_ex_run_coloc` | GPL | Colocalization analysis for shared genetic associations |
| **hyprcoloc** | `util_ex_run_hyprcoloc` | GPL-2.0-or-later | Hierarchical Bayesian colocalization analysis |
| **TwoSampleMR** | `util_ex_run_2samplemr` | MIT | Mendelian Randomization using two-sample summary data |
| **MESuSiE** | `util_ex_run_mesusie` | GPL (>= 2) | Multivariate extension of SuSiE for fine-mapping. R package by borangao. Original reference: Gao & Zhou, Nature Genetics, 2024. |
| **PLINK/PLINK2** | `util_ex_run_clumping`, `util_ex_calculate_ldmatrix`, `util_ex_process_ref` | GPL-3.0 | PLINK whole genome association analysis toolset |
| **BCFtools** | `hm_assign_rsid`, `util_ex_process_h5` | MIT/Expat License | Variant calling and manipulation utilities. Part of samtools/htslib |
| **MAGMA** | `util_ex_run_magma` | GPL v3.0 (v1.0 only); Proprietary (v1.01+, free for academic use) | Multi-marker Analysis of GenoMic Annotation. Original v1.0 is GPL v3.0; subsequent versions use standard copyright |
| **tabix** | `util_ex_ldproxyfinder` | MIT License | Part of samtools/htslib |
| **SCDRS** | `util_ex_run_scdrs` | MIT | Single-cell Deconvolution and Regression for Summary statistics |

**Note on R Packages**: 
- **susieR** uses a BSD-style license (permissive, GPL-3.0 compatible)
- **coloc** uses GPL license (GPL-3.0 compatible)
- **hyprcoloc** is licensed under GPL-2.0-or-later (GPL-3.0 compatible)
- **MESuSiE** is licensed under GPL (>= 2), equivalent to GPL-2.0-or-later (GPL-3.0 compatible)
- **TwoSampleMR** uses MIT license (permissive, GPL-3.0 compatible)
- These packages are called via R subprocess calls and are not included in the GWASLab codebase

**Important**: These external tools have their own licenses and installation requirements. Users must install and license these tools separately if they wish to use the corresponding GWASLab functions. GWASLab's GPL-3.0 license does not apply to these external tools.

## License Compatibility Analysis

### Summary

All dependencies use permissive licenses (BSD-3-Clause, MIT, Apache-2.0, PSF) that are compatible with GPL-3.0. This means:

1. **GWASLab can remain GPL-3.0-only**: All dependencies are compatible with GPL-3.0, so the project can maintain its GPL-3.0-only license.

2. **No license conflicts**: There are no copyleft dependencies that would require GWASLab to change its license.

3. **Distribution compliance**: When distributing GWASLab, you must:
   - Include the GPL-3.0 license text (already in root `LICENSE` file)
   - Provide source code or make it available
   - Include license notices for all dependencies (this document serves that purpose)

### License Types Used

- **GPL-3.0-only**: GWASLab itself (copyleft license)
- **BSD-3-Clause**: pandas, numpy, seaborn, scipy, h5py (permissive)
- **BSD-style**: susieR (permissive, GPL-compatible)
- **MIT**: pysam, adjustText, scikit-allel, polars, sumstats-liftover, setuptools, wheel, TwoSampleMR (permissive)
- **Apache-2.0**: pyarrow, requests (permissive)
- **PSF**: matplotlib (permissive, GPL-compatible)
- **GPL**: coloc (copyleft, GPL-3.0 compatible)
- **GPL-2.0-or-later**: hyprcoloc (copyleft, GPL-3.0 compatible)
- **GPL (>= 2)**: MESuSiE (copyleft, GPL-3.0 compatible, equivalent to GPL-2.0-or-later)
- **GPL-3.0**: MultiSuSiE, LDSC (copyleft, GPL-3.0 compatible)
- **GPL v3.0 (v1.0 only) / Proprietary**: MAGMA (v1.0 is GPL v3.0; later versions are proprietary)

### Important Notes

1. **GPL-3.0 Copyleft**: Since GWASLab is GPL-3.0-only, any derivative works or modifications must also be licensed under GPL-3.0 or a compatible license.

2. **Dependency Licenses**: The permissive licenses (BSD, MIT, Apache-2.0, PSF) allow their code to be used in GPL-licensed projects. However, the GPL-3.0 license of GWASLab applies to the combined work.

3. **Distribution Requirements**: When distributing GWASLab:
   - The full GPL-3.0 license must be included
   - Source code must be made available
   - All copyright notices from dependencies must be preserved
   - This license summary should be included

## References

- [GPL-3.0 License](https://www.gnu.org/licenses/gpl-3.0.html)
- [BSD-3-Clause License](https://opensource.org/licenses/BSD-3-Clause)
- [MIT License](https://opensource.org/licenses/MIT)
- [Apache-2.0 License](https://www.apache.org/licenses/LICENSE-2.0)
- [Python Software Foundation License](https://docs.python.org/3/license.html)

## Last Updated

This document was generated by analyzing the project dependencies listed in `pyproject.toml` and `environment.yml`. License information was retrieved from PyPI and package repositories.

For the most up-to-date license information, please refer to:
- Each package's official repository
- PyPI package pages: https://pypi.org/
- Package LICENSE files in their source distributions

