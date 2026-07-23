#!/usr/bin/env bash
# =============================================================================
# 10_fill.sh — Statistics fill (--fill)
# =============================================================================
# Derive missing statistical columns via fill_data() on the CLI.
# See docs/Conversion.md for conversion formulas and priority.
#
# Usage:
#   bash 10_fill.sh
#
# Requires:
#   pip install gwaslab
# =============================================================================

set -euo pipefail

mkdir -p output

# Minimal MLOG10P-only fixture for demonstration
ML10P_ONLY="output/10_mlog10p_only.tsv"
cat > "$ML10P_ONLY" <<'EOF'
SNPID	CHR	POS	EA	NEA	MLOG10P
rs1	1	100	A	G	8.0
rs2	1	200	G	A	10.0
rs3	2	300	C	T	5.5
EOF

echo "=== [1] MLOG10P → P, export gwaslab format ==="
gwaslab \
    --input   "$ML10P_ONLY" \
    --fill    P \
    --no-gzip \
    --output  "output/10_with_p.gwaslab.tsv" \
    --to-fmt  gwaslab

echo ""
echo "=== [2] MLOG10P → P, export LDSC format ==="
gwaslab \
    --input   "$ML10P_ONLY" \
    --fill    P \
    --no-gzip \
    --output  "output/10_with_p.ldsc" \
    --to-fmt  ldsc

echo ""
echo "=== [3] EAF → MAF, then MAF filter ==="
gwaslab \
    --input   "../../test/raw/dirty_sumstats.tsv" \
    --fill    MAF \
    --maf     0.01 \
    --no-gzip \
    --output  "output/10_maf_filtered.tsv" \
    --to-fmt  gwaslab

echo ""
echo "Done. Output files written under output/"
