#!/usr/bin/env python3
"""Print runtime verification for good_patterns.py (run after compileall check)."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import good_patterns as gp  # noqa: E402


def main() -> int:
    checks: list[tuple[str, bool]] = []

    checks.append(("sep value", gp.SEP_GOOD == gp.EXPECTED_SEP))
    checks.append(("ylabel value", gp.YLABEL_GOOD == gp.EXPECTED_YLABEL))
    checks.append(
        (
            "ldsc regex",
            gp._LDSC_LINE_RE.findall(gp.SAMPLE_LDSC_LINE) == ["obs:", "0.25", "se:", "0.03", "NA"],
        )
    )
    checks.append(("status regex", gp.STATUS_RE.match(gp.SAMPLE_STATUS) is not None))
    checks.append(
        (
            "legend dynamic",
            gp.legend_title_good(5.0, -8) == "$\\mathregular{ P < 5.0 x 10^{-8} }$ in:",
        )
    )
    checks.append(
        (
            "p text dynamic",
            gp.p_text_good("1.2", "-8") == "$\\mathregular{p = 1.2 \\times  10^{-8}}$",
        )
    )
    checks.append(
        (
            "plink script newlines",
            "\\\n" in gp.plink_script_good("--bfile x", "out")
            and "plink2" in gp.plink_script_good("--bfile x", "out"),
        )
    )
    checks.append(
        (
            "docstring roundtrip",
            gp.parse_pattern_doc_example(r"\w\w\w\w\w[35]\w") == r"\w\w\w\w\w[35]\w",
        )
    )

    failed = [name for name, ok in checks if not ok]
    for name, ok in checks:
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}")

    if failed:
        print(f"\n{len(failed)} check(s) failed.", file=sys.stderr)
        return 1

    print(f"\nAll {len(checks)} checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
