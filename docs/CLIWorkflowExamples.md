# CLI Workflow Examples

This document provides complete workflow examples for using the GWASLab command-line interface. For basic usage and argument reference, see [CLI](CLI.md).

## Example 1: Basic QC Pipeline

```python
# Load, QC, and save sumstats
gwaslab --input raw_sumstats.tsv \
  --fmt auto \
  --qc \
  --remove \
  --remove-dup \
  --normalize \
  --threads 4 \
  --out cleaned_sumstats \
  --to-fmt gwaslab \
  --quiet
```

This will:
1. Load sumstats from `raw_sumstats.tsv` with auto-detection
2. Perform QC (fix IDs, chromosomes, positions, alleles)
3. Remove bad quality variants
4. Remove duplicates
5. Normalize indels
6. Save to `cleaned_sumstats.gwaslab.tsv.gz`

## Example 2: Harmonization Pipeline

```python
# Full harmonization with all reference files
gwaslab --input sumstats.tsv \
  --fmt auto \
  --harmonize \
  --ref-seq /data/reference/hg19.fasta.gz \
  --ref-rsid-vcf /data/reference/dbsnp.vcf.gz \
  --ref-infer /data/reference/1kg.vcf.gz \
  --ref-alt-freq AF \
  --maf-threshold 0.40 \
  --ref-maf-threshold 0.5 \
  --sweep-mode \
  --threads 8 \
  --out harmonized_sumstats \
  --to-fmt gwaslab \
  --tab-fmt tsv
```

This will:
1. Load sumstats
2. Run basic QC
3. Check alignment with reference sequence
4. Flip alleles if needed
5. Assign rsIDs from VCF
6. Infer strand for palindromic SNPs
7. Save harmonized sumstats

## Example 3: Format Conversion for LDSC

```python
# Convert to LDSC format for heritability analysis
gwaslab --input sumstats.tsv \
  --fmt gwaslab \
  --out sumstats_ldsc \
  --to-fmt ldsc \
  --no-gzip
```

## Example 4: Prepare Data for Meta-analysis

```python
# QC, harmonize, and format for meta-analysis
gwaslab --input study1.tsv \
  --fmt auto \
  --qc \
  --remove-dup \
  --harmonize \
  --ref-seq /data/reference/hg19.fasta.gz \
  --ref-rsid-vcf /data/reference/dbsnp.vcf.gz \
  --out study1_ready \
  --to-fmt gwaslab \
  --exclude-hla \
  --n 50000
```

## Example 5: Extract HapMap3 Variants

```python
# Extract only HapMap3 variants for replication
gwaslab --input discovery_sumstats.tsv \
  --fmt gwaslab \
  --out hapmap3_sumstats \
  --to-fmt gwaslab \
  --hapmap3 \
  --build 19
```

## Example 6: Multi-step Processing

```python
# Step 1: QC
gwaslab --input raw.tsv --fmt auto --qc --remove-dup --out step1_qc --to-fmt gwaslab

# Step 2: Harmonize
gwaslab --input step1_qc.gwaslab.tsv.gz --fmt gwaslab --harmonize \
  --ref-seq ref.fasta --out step2_harmonized --to-fmt gwaslab

# Step 3: Format for specific tool
gwaslab --input step2_harmonized.gwaslab.tsv.gz --fmt gwaslab \
  --out final --to-fmt plink --tab-fmt tsv
```

## Example 7: Assign rsID Only

```python
# Assign rsID from VCF file
gwaslab --input sumstats.tsv \
  --fmt auto \
  --assign-rsid \
  --ref-rsid-vcf /path/to/dbsnp.vcf.gz \
  --overwrite empty \
  --threads 4 \
  --out sumstats_with_rsid \
  --to-fmt gwaslab
```

## Example 8: Convert rsID to CHR:POS

```python
# Convert rsID to CHR:POS using VCF (auto-generates HDF5)
gwaslab --input sumstats.tsv \
  --fmt auto \
  --rsid-to-chrpos \
  --ref-rsid-vcf /path/to/reference.vcf.gz \
  --build 19 \
  --threads 8 \
  --out sumstats_with_chrpos \
  --to-fmt gwaslab
```

## Example 9: Complete Pipeline with All Options

```python
# Complete pipeline: QC + Harmonization + Format conversion
gwaslab --input raw_sumstats.tsv \
  --fmt auto \
  --qc \
  --remove \
  --remove-dup \
  --normalize \
  --harmonize \
  --ref-seq /data/reference/hg19.fasta.gz \
  --ref-rsid-vcf /data/reference/dbsnp.vcf.gz \
  --ref-infer /data/reference/1kg.vcf.gz \
  --ref-alt-freq AF \
  --maf-threshold 0.40 \
  --ref-maf-threshold 0.5 \
  --sweep-mode \
  --threads 8 \
  --out final_sumstats \
  --to-fmt ldsc \
  --tab-fmt tsv \
  --exclude-hla \
  --hapmap3 \
  --n 100000 \
  --chr-prefix chr \
  --xymt-number
```

## Example 10: Testing with Small Sample

```python
# Test command on first 1000 rows
gwaslab --input large_sumstats.tsv \
  --fmt auto \
  --nrows 1000 \
  --qc \
  --remove-dup \
  --out test_output \
  --to-fmt gwaslab
```

## Example 11: Output with bgzip and tabix

```python
# Create bgzip-compressed file with tabix index
gwaslab --input sumstats.tsv \
  --fmt gwaslab \
  --out indexed_sumstats \
  --to-fmt gwaslab \
  --bgzip \
  --tabix
```

## Example 12: Custom HLA Region Exclusion

```python
# Exclude custom HLA region (chr6:20-30 Mbp)
gwaslab --input sumstats.tsv \
  --fmt gwaslab \
  --out no_hla \
  --to-fmt gwaslab \
  --exclude-hla \
  --hla-lower 20 \
  --hla-upper 30
```

