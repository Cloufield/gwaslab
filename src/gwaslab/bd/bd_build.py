"""
Genome build mappings and display names for different species.

This module provides species-specific genome build mappings that normalize
various build identifier formats to standardized two-digit numeric codes.
All build codes have a non-zero first digit (e.g., "19", "20", "36").
"""

# Species-specific genome build mappings
# Format: {species: {build_aliases: two_digit_build_code, ...}}
# All builds are mapped to two-digit numeric strings where the first digit is never 0
# (e.g., "19", "20", "36", "46" - not "09", "06", etc.)
BUILD_MAPPINGS = {
    "homo sapiens": {
        ("hg19", "19", "37", "b37", "grch37", "grch37.p13"): "19",
        ("hg18", "18", "36", "b36", "grch36"): "18",
        ("hg38", "38", "b38", "grch38", "grch38.p14"): "38",
        ("t2t", "hs1", "chm13", "13", "t2t-chm13"): "13",
    },
    "human": {
        ("hg19", "19", "37", "b37", "grch37", "grch37.p13"): "19",
        ("hg18", "18", "36", "b36", "grch36"): "18",
        ("hg38", "38", "b38", "grch38", "grch38.p14"): "38",
        ("t2t", "hs1", "chm13", "13", "t2t-chm13"): "13",
    },
    "mus musculus": {
        ("mm10", "10", "grcm38"): "20",
        ("mm9", "9", "ncbi37"): "29",
        ("mm39", "39", "grcm39"): "39",
    },
    "mouse": {
        ("mm10", "10", "grcm38"): "20",
        ("mm9", "9", "ncbi37"): "29",
        ("mm39", "39", "grcm39"): "39",
    },
    "rattus norvegicus": {
        ("rn6", "6", "rgsc6.0"): "36",
        ("rn7", "7", "mratbn7.2"): "37",
        ("rn5", "5", "rgsc5.0"): "35",
    },
    "rat": {
        ("rn6", "6", "rgsc6.0"): "36",
        ("rn7", "7", "mratbn7.2"): "37",
        ("rn5", "5", "rgsc5.0"): "35",
    },
    "gallus gallus": {
        ("galgal6", "galgal6a", "6", "grcg6a"): "46",
        ("galgal5", "5", "grcg5"): "45",
        ("galgal4", "4"): "44",
    },
    "chicken": {
        ("galgal6", "galgal6a", "6", "grcg6a"): "46",
        ("galgal5", "5", "grcg5"): "45",
        ("galgal4", "4"): "44",
    },
    "danio rerio": {
        ("zv11", "11", "grcz11"): "51",
        ("zv10", "10", "grcz10"): "50",
        ("zv9", "9", "grcz9"): "59",
    },
    "zebrafish": {
        ("zv11", "11", "grcz11"): "51",
        ("zv10", "10", "grcz10"): "50",
        ("zv9", "9", "grcz9"): "59",
    },
    "sus scrofa": {
        ("susscr11", "11", "ssc11"): "61",
        ("susscr10", "10", "ssc10"): "60",
    },
    "pig": {
        ("susscr11", "11", "ssc11"): "61",
        ("susscr10", "10", "ssc10"): "60",
    },
    "bos taurus": {
        ("ars-ucd1.2", "ars1.2", "ars_ucd1.2"): "72",
        ("umd3.1", "3.1"): "73",
    },
    "cattle": {
        ("ars-ucd1.2", "ars1.2", "ars_ucd1.2"): "72",
        ("umd3.1", "3.1"): "73",
    },
    "cow": {
        ("ars-ucd1.2", "ars1.2", "ars_ucd1.2"): "72",
        ("umd3.1", "3.1"): "73",
    },
    "canis lupus familiaris": {
        ("canfam3.1", "3.1"): "83",
        ("canfam2", "2"): "82",
    },
    "dog": {
        ("canfam3.1", "3.1"): "83",
        ("canfam2", "2"): "82",
    },
    "equus caballus": {
        ("equcab3.0", "3.0", "ecab3"): "93",
        ("equcab2.0", "2.0"): "92",
    },
    "horse": {
        ("equcab3.0", "3.0", "ecab3"): "93",
        ("equcab2.0", "2.0"): "92",
    },
    "oryza sativa": {
        ("irgsp-1.0", "irgsp1.0", "1.0"): "81",
        ("msu7", "7"): "87",
    },
    "rice": {
        ("irgsp-1.0", "irgsp1.0", "1.0"): "81",
        ("msu7", "7"): "87",
    },
    "arabidopsis thaliana": {
        ("tair10", "10"): "91",
        ("tair9", "9"): "99",
    },
    "arabidopsis": {
        ("tair10", "10"): "91",
        ("tair9", "9"): "99",
    },
}

# Build display names for logging (species-specific mappings)
# Format: {(species, build_code): display_name}
BUILD_DISPLAY_NAMES = {
    # Human builds
    ("homo sapiens", "19"): "GRCh37/hg19",
    ("human", "19"): "GRCh37/hg19",
    ("homo sapiens", "18"): "GRCh36/hg18",
    ("human", "18"): "GRCh36/hg18",
    ("homo sapiens", "38"): "GRCh38/hg38",
    ("human", "38"): "GRCh38/hg38",
    ("homo sapiens", "13"): "T2T-CHM13",
    ("human", "13"): "T2T-CHM13",
    # Mouse builds
    ("mus musculus", "20"): "GRCm38/mm10",
    ("mouse", "20"): "GRCm38/mm10",
    ("mus musculus", "29"): "NCBI37/mm9",
    ("mouse", "29"): "NCBI37/mm9",
    ("mus musculus", "39"): "GRCm39/mm39",
    ("mouse", "39"): "GRCm39/mm39",
    # Rat builds
    ("rattus norvegicus", "36"): "RGSC6.0/rn6",
    ("rat", "36"): "RGSC6.0/rn6",
    ("rattus norvegicus", "37"): "mRatBN7.2/rn7",
    ("rat", "37"): "mRatBN7.2/rn7",
    ("rattus norvegicus", "35"): "RGSC5.0/rn5",
    ("rat", "35"): "RGSC5.0/rn5",
    # Chicken builds
    ("gallus gallus", "46"): "GRCg6a/galGal6",
    ("chicken", "46"): "GRCg6a/galGal6",
    ("gallus gallus", "45"): "GRCg5/galGal5",
    ("chicken", "45"): "GRCg5/galGal5",
    ("gallus gallus", "44"): "galGal4",
    ("chicken", "44"): "galGal4",
    # Zebrafish builds
    ("danio rerio", "51"): "GRCz11/zv11",
    ("zebrafish", "51"): "GRCz11/zv11",
    ("danio rerio", "50"): "GRCz10/zv10",
    ("zebrafish", "50"): "GRCz10/zv10",
    ("danio rerio", "59"): "GRCz9/zv9",
    ("zebrafish", "59"): "GRCz9/zv9",
    # Pig builds
    ("sus scrofa", "61"): "Sscrofa11.1/susScr11",
    ("pig", "61"): "Sscrofa11.1/susScr11",
    ("sus scrofa", "60"): "Sscrofa10.2/susScr10",
    ("pig", "60"): "Sscrofa10.2/susScr10",
    # Cattle builds
    ("bos taurus", "72"): "ARS-UCD1.2",
    ("cattle", "72"): "ARS-UCD1.2",
    ("cow", "72"): "ARS-UCD1.2",
    ("bos taurus", "73"): "UMD3.1",
    ("cattle", "73"): "UMD3.1",
    ("cow", "73"): "UMD3.1",
    # Dog builds
    ("canis lupus familiaris", "83"): "CanFam3.1",
    ("dog", "83"): "CanFam3.1",
    ("canis lupus familiaris", "82"): "CanFam2",
    ("dog", "82"): "CanFam2",
    # Horse builds
    ("equus caballus", "93"): "EquCab3.0",
    ("horse", "93"): "EquCab3.0",
    ("equus caballus", "92"): "EquCab2.0",
    ("horse", "92"): "EquCab2.0",
    # Rice builds
    ("oryza sativa", "81"): "IRGSP-1.0",
    ("rice", "81"): "IRGSP-1.0",
    ("oryza sativa", "87"): "MSU7",
    ("rice", "87"): "MSU7",
    # Arabidopsis builds
    ("arabidopsis thaliana", "91"): "TAIR10",
    ("arabidopsis", "91"): "TAIR10",
    ("arabidopsis thaliana", "99"): "TAIR9",
    ("arabidopsis", "99"): "TAIR9",
}

