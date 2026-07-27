"""Fixed escape-sequence patterns — must compile with -W error::SyntaxWarning."""

from __future__ import annotations

import re

# --- Category A: regex / sep (ldsc_parse, io_read_ldsc, qc_fix_sumstats_polars, ...) ---

SEP_GOOD = r"\s+"
_LDSC_LINE_RE = re.compile(r"[a-zA-Z\s\d]+:|[-0-9.]+[e]?[-0-9.]+|NA")
_SFILE_RE = re.compile(r"--sfile1 ([^\s]+) --sfile2 ([^\s]+)[ /n]")
STATUS_RE = re.compile(r"^\w\w\w\w\w[45]\w")
REPLACE_GOOD = "file@name".replace("@", r"(\w+)")

# --- Category B: matplotlib LaTeX (viz_plot_miamiplot2, viz_plot_compare_effect, ...) ---

YLABEL_GOOD = r"$\mathregular{-log_{10}(P)}$"


def legend_title_good(mantissa: float, exponent: int) -> str:
    return rf"$\mathregular{{ P < {mantissa} x 10^{{{exponent}}} }}$ in:"


def p_text_good(p12: str, pe: str) -> str:
    return rf"$\mathregular{{p = {p12} \times  10^{{{pe}}}}}$"


# --- Category C: docstring examples (g_vchange_status, io_load_ld, multisusie) ---


def parse_pattern_doc_example(pattern: str) -> str:
    r"""Document regex pattern format.

    Parameters
    ----------
    pattern : str
        Regex pattern like r'\w\w\w\w\w[35]\w'

    Returns
    -------
    str
        Echo of *pattern*.
    """
    return pattern


DOC_GOOD = r"""
Default separator is whitespace (``\s+``).
LaTeX example: $\\sigma^2_k$
Regex example: r'\w\w\w\w\w[35]\w'
"""


# --- Category D: shell script (util_ex_plink_filter) ---


def plink_script_good(filter_flag: str, out_prefix: str) -> str:
    return (
        "plink2 \\\n"
        f"    {filter_flag} \\\n"
        "    --make-just-bim \\\n"
        "    --make-just-fam \\\n"
        f"    --out {out_prefix}\n"
    )


# --- Expected runtime values (same as deprecated bad literals) ---

EXPECTED_SEP = "\\s+"
EXPECTED_YLABEL = "$\\mathregular{-log_{10}(P)}$"
SAMPLE_LDSC_LINE = "h2_obs: 0.25 h2_se: 0.03 NA"
SAMPLE_STATUS = "1234547"
