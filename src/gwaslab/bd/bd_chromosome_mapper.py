"""
Unified ChromosomeMapper class for handling chromosome notation conversions.

This module provides a single, comprehensive interface for mapping chromosome
identifiers between sumstats and reference files using a numeric middle layer.
"""

import pandas as pd
import numpy as np
from typing import Union, Optional, Dict, Any, Tuple
import re
import os

from gwaslab.bd.bd_sex_chromosomes import Chromosomes
from gwaslab.bd.bd_common_data import get_chr_to_NC
from gwaslab.info.g_Log import Log


class ChromosomeMapper:
    """
    ChromosomeMapper with numeric middle layer architecture (species-aware).
    
    Architecture:
    - Middle layer: Numbers (species-specific, e.g., 1-22, 23, 24, 25 for human X, Y, MT)
    - Sumstats layer: Maps sumstats chr notation → numbers → sumstats chr notation
    - Reference layer: Maps reference chr notation → numbers → reference chr notation
    
    This design makes it easy to match sumstats chr with reference chr:
    - sumstats_chr → number → reference_chr
    - reference_chr → number → sumstats_chr
    
    The mapper is fully species-aware and handles:
    - Different chromosome counts per species (e.g., human: 22 autosomes, mouse: 19)
    - Different sex chromosome notation (e.g., X/Y for mammals, Z/W for birds)
    - Species-specific mitochondrial notation
    
    **Unconvertible Chromosomes:**
    
    When a chromosome identifier cannot be converted to a number, the mapper returns
    `pd.NA` instead of raising an error. This includes:
    - Alternative loci (e.g., "1_KI270766v1_alt", "KI270766.1")
    - Unplaced sequences (e.g., "GL000195.1", "GL000009.2", "chrUn_GL000195v1")
    - Chromosome arm notation (e.g., "1p", "1q", "chr1p")
    - Malformed or invalid identifiers
    - Out-of-range chromosome numbers
    
    This allows processing of mixed datasets with both standard and non-standard
    chromosome identifiers without interrupting the workflow.
    
    Parameters
    ----------
    species : str, default="homo sapiens"
        Species name (case-insensitive). Uses Chromosomes class for species-specific mappings.
        Supported species include: homo sapiens, mus musculus, rattus norvegicus,
        gallus gallus, danio rerio, drosophila melanogaster, sus scrofa, bos taurus,
        canis lupus familiaris, equus caballus, oryza sativa, arabidopsis thaliana.
    build : str, optional
        Genome build version (e.g., "19", "38") for NCBI RefSeq ID mappings.
        If not provided or unavailable for the species/build combination, NCBI RefSeq
        conversions will not be available (no warnings are logged).
        
        **Build Auto-Detection:**
        - For **reference files** with NCBI RefSeq notation (NC_*), the build can be
          automatically detected by matching contig IDs against known builds (hg19/hg38).
          This happens when `detect_reference_format()` is called on a file containing
          NC notation and build is not set.
        - For **sumstats data** with NCBI RefSeq notation, the build parameter **must**
          be explicitly specified. Build is not auto-detected from sumstats data alone.
    xymt_num : list, default=[23, 24, 25]
        Numeric values for sex chromosomes and mitochondrial (human convention).
        For other species, adjust based on their chromosome structure.
    log : Log, optional
        Logging object for warnings and messages.
    verbose : bool, default=True
        Whether to print warnings and messages.
    
    Notes
    -----
    **Behavior without `detect_sumstats_format()`:**
    
    The mapper defaults to **numeric format** if `detect_sumstats_format()` has not
    been called. This means:
    
    - `sumstats_to_number()` and `number_to_sumstats()` will automatically build
      a numeric sumstats layer on first use
    - Numeric chromosomes (1, 2, 23) work directly
    - String chromosomes ("X", "Y", "MT") use fallback string mappings
    - Format-specific strings like "chr1" may not work correctly without explicit
      format detection
    
    It is **recommended** to call `detect_sumstats_format()` to ensure all
    format-specific conversions work correctly for your data format.
    
    **Return Values:**
    
    Methods that convert chromosomes to numbers (e.g., `sumstats_to_number()`,
    `reference_to_number()`, `to_numeric()`) return:
    - `int` for successfully converted standard chromosomes
    - `pd.NA` for unconvertible chromosomes (alternative loci, unplaced sequences,
      malformed identifiers, etc.)
    
    When working with pandas Series, unconvertible values will appear as `pd.NA`
    in the resulting Series, allowing you to filter or handle them as needed.
    
    **Supported Reference File Formats:**
    
    The mapper can automatically detect chromosome format from:
    - VCF files (`.vcf`, `.vcf.gz`, `.bcf`)
    - FASTA files (`.fa`, `.fasta`, `.fa.gz`, `.fasta.gz`, `.fa.bgz`, `.fasta.bgz`)
    - GTF/GFF files (`.gtf`, `.gff`, `.gtf.gz`, `.gff.gz`)
    - Chain files (`.chain`, `.chain.gz`)
    
    Format detection is case-insensitive and handles various chromosome notations:
    - Numeric: "1", "2", "22", "23", "24", "25"
    - String: "1", "2", "X", "Y", "MT"
    - Chr-prefixed: "chr1", "chr2", "chrX", "chrY", "chrMT" (case variations supported)
    - NCBI RefSeq: "NC_000001.11", "NC_000023.11"
      * For **reference files**: Build can be auto-detected from NC notation
      * For **sumstats data**: Build parameter must be explicitly specified
      * For **reference files**: Build can be auto-detected from NC notation
      * For **sumstats data**: Build parameter must be explicitly specified
    
    Examples
    --------
    Basic usage with auto-detection:
    
    >>> import pandas as pd
    >>> from gwaslab.bd.bd_chromosome_mapper import ChromosomeMapper
    >>> 
    >>> # Initialize mapper for human (default)
    >>> mapper = ChromosomeMapper(species="homo sapiens")
    >>> 
    >>> # Auto-detect sumstats format and build sumstats layer
    >>> sumstats_chr = pd.Series([1, 2, 23, "X", "MT"])
    >>> mapper.detect_sumstats_format(sumstats_chr)
    'numeric'
    >>> 
    >>> # Convert sumstats chromosomes to reference format (auto-detects from file)
    >>> mapper.sumstats_to_reference(1, reference_file="reference.vcf.gz")
    'chr1'
    >>> mapper.sumstats_to_reference("X", reference_file="reference.vcf.gz")
    'chrX'
    >>> 
    >>> # Convert entire Series
    >>> mapper.sumstats_to_reference_series(sumstats_chr, reference_file="reference.vcf.gz")
    0       chr1
    1       chr2
    2      chr23
    3       chrX
    4      chrMT
    dtype: object
    
    Handling unconvertible chromosomes:
    
    >>> # Unconvertible chromosomes return pd.NA
    >>> mapper.sumstats_to_number("1_KI270766v1_alt")
    <NA>
    >>> mapper.sumstats_to_number("GL000195.1")
    <NA>
    >>> mapper.sumstats_to_number("1p")  # Chromosome arm notation
    <NA>
    >>> 
    >>> # Series with mixed valid and invalid chromosomes
    >>> series = pd.Series(["1", "2", "1_KI270766v1_alt", "X", "GL000195.1"])
    >>> result = mapper.to_numeric(series)
    >>> result
    0       1
    1       2
    2    <NA>
    3      23
    4    <NA>
    dtype: object
    
    Working with different sumstats formats:
    
    >>> # Numeric format (1, 2, 23, 24, 25)
    >>> mapper = ChromosomeMapper()
    >>> mapper.detect_sumstats_format(pd.Series([1, 2, 23]))
    'numeric'
    >>> mapper.sumstats_to_number(1)
    1
    >>> mapper.number_to_sumstats(23)
    23
    
    >>> # String format ("1", "2", "X", "Y", "MT")
    >>> mapper.detect_sumstats_format(pd.Series(["1", "2", "X", "MT"]))
    'string'
    >>> mapper.sumstats_to_number("X")
    23
    >>> mapper.number_to_sumstats(23)
    'X'
    
    >>> # Chr-prefixed format ("chr1", "chr2", "chrX")
    >>> mapper.detect_sumstats_format(pd.Series(["chr1", "chr2", "chrX"]))
    'chr'
    >>> mapper.sumstats_to_number("chr1")
    1
    >>> mapper.number_to_sumstats(1)
    'chr1'
    
    >>> # NCBI RefSeq format in sumstats (requires build parameter)
    >>> mapper = ChromosomeMapper(build="38")
    >>> mapper.detect_sumstats_format(pd.Series(["NC_000001.11", "NC_000023.11"]))
    'nc'
    >>> mapper.sumstats_to_number("NC_000001.11")
    1
    >>> mapper.number_to_sumstats(1)
    'NC_000001.11'
    
    >>> # NCBI RefSeq format in reference file (build auto-detected)
    >>> mapper = ChromosomeMapper()  # build not specified
    >>> mapper.detect_reference_format("reference.vcf.gz")  # Contains NC_000001.11
    ('nc', '')
    >>> # Build is auto-detected as "38" from NC notation in the file
    >>> mapper.sumstats_to_reference(1, reference_file="reference.vcf.gz")
    'NC_000001.11'
    
    Matching sumstats with reference files:
    
    >>> # Sumstats uses numeric format, reference VCF uses "chr1", "chr2"
    >>> mapper = ChromosomeMapper()
    >>> mapper.detect_sumstats_format(pd.Series([1, 2, 23]))
    >>> mapper.sumstats_to_reference(1, reference_file="reference.vcf.gz")
    'chr1'
    >>> mapper.sumstats_to_reference(23, reference_file="reference.vcf.gz")
    'chrX'
    >>> 
    >>> # Reverse: reference uses "chr1", convert to sumstats format
    >>> mapper.reference_to_sumstats("chr1", reference_file="reference.vcf.gz")
    1
    >>> mapper.reference_to_sumstats("chrX", reference_file="reference.vcf.gz")
    23
    
    Using convenience methods:
    
    >>> mapper = ChromosomeMapper()
    >>> mapper.detect_sumstats_format(pd.Series(["chr1", "chr2", "chrX"]))
    >>> 
    >>> # Convert to numeric (middle layer)
    >>> mapper.to_numeric("chr1")
    1
    >>> mapper.to_numeric(pd.Series(["chr1", "chr2", "chrX"]))
    0     1
    1     2
    2    23
    dtype: int64
    >>> 
    >>> # Convert to string format
    >>> mapper.to_string(1)
    '1'
    >>> mapper.to_string(23)
    'X'
    >>> 
    >>> # Convert to chr-prefixed format
    >>> mapper.to_chr(1)
    'chr1'
    >>> mapper.to_chr(23)
    'chrX'
    >>> 
    >>> # Convert to NCBI RefSeq (requires build)
    >>> mapper = ChromosomeMapper(build="38")
    >>> mapper.detect_sumstats_format(pd.Series([1, 2, 23]))
    >>> mapper.to_nc(1)
    'NC_000001.11'
    >>> mapper.to_nc(23)
    'NC_000023.11'
    
    Different species support:
    
    >>> # Human (default): 22 autosomes, X, Y, MT
    >>> mapper = ChromosomeMapper(species="homo sapiens")
    >>> mapper.detect_sumstats_format(pd.Series([1, 2, 22, 23, 24, 25]))
    'numeric'
    >>> 
    >>> # Mouse: 19 autosomes, X, Y, MT
    >>> mapper_mouse = ChromosomeMapper(species="mus musculus")
    >>> mapper_mouse.detect_sumstats_format(pd.Series([1, 2, 19, 20, 21, 22]))
    'numeric'
    >>> 
    >>> # Chicken: uses Z/W instead of X/Y
    >>> mapper_chicken = ChromosomeMapper(species="gallus gallus")
    >>> mapper_chicken.detect_sumstats_format(pd.Series(["1", "2", "Z", "W", "MT"]))
    'string'
    >>> mapper_chicken.sumstats_to_number("Z")
    20
    >>> mapper_chicken.number_to_sumstats(20)
    'Z'
    
    Integration with Sumstats object:
    
    >>> from gwaslab.g_Sumstats import Sumstats
    >>> 
    >>> # Sumstats object automatically creates and manages a mapper
    >>> mysumstats = Sumstats(data, verbose=False)
    >>> 
    >>> # Mapper auto-detects format from data
    >>> print(mysumstats.mapper._sumstats_format)
    'numeric'
    >>> 
    >>> # Use mapper for reference file matching
    >>> mysumstats.check_ref(ref_seq="reference.fa", verbose=False)
    >>> # Mapper automatically handles chromosome conversion between sumstats and FASTA
    
    Working with reference files (auto-detection):
    
    >>> mapper = ChromosomeMapper()
    >>> mapper.detect_sumstats_format(pd.Series([1, 2, 23]))
    >>> 
    >>> # VCF file with "chr1", "chr2" format
    >>> mapper.detect_reference_format("reference.vcf.gz")
    >>> mapper.sumstats_to_reference(1, reference_file="reference.vcf.gz")
    'chr1'
    >>> 
    >>> # FASTA file with "1", "2" format
    >>> mapper.detect_reference_format("reference.fa")
    >>> mapper.sumstats_to_reference(1, reference_file="reference.fa")
    1
    >>> 
    >>> # GTF file with "chr1", "chr2" format
    >>> mapper.detect_reference_format("annotation.gtf.gz")
    >>> mapper.sumstats_to_reference(1, reference_file="annotation.gtf.gz")
    'chr1'
    
    >>> # Reference file with NCBI RefSeq notation (build auto-detected)
    >>> mapper = ChromosomeMapper()  # build not specified
    >>> mapper.detect_reference_format("reference_nc.vcf.gz")  # Contains NC_000001.11
    ('nc', '')
    >>> # Build "38" is automatically detected from NC notation
    >>> mapper.sumstats_to_reference(1, reference_file="reference_nc.vcf.gz")
    'NC_000001.11'
    """
    
    # Format detection patterns
    _NCBI_PATTERN = re.compile(r'^(NC_|CM_)\d+\.\d+$', re.IGNORECASE)
    _CHR_PREFIX_PATTERN = re.compile(r'^(chr|Chr|CHR)', re.IGNORECASE)
    
    # ============================================================================
    # Initialization
    # ============================================================================
    
    def __init__(self,
                 species: str = "homo sapiens",
                 build: Optional[str] = None,
                 xymt_num: list = [23, 24, 25],
                 log: Optional[Log] = None,
                 verbose: bool = True):
        """
        Initialize ChromosomeMapper with numeric middle layer architecture.
        
        Parameters
        ----------
        species : str, default="homo sapiens"
            Species name (case-insensitive). The mapper will automatically use
            species-specific chromosome counts and notation (e.g., X/Y for mammals,
            Z/W for birds). See class docstring for supported species.
        build : str, optional
            Genome build version (e.g., "19", "38") for NCBI RefSeq ID mappings.
        xymt_num : list, default=[23, 24, 25]
            Numeric values for sex chromosomes and mitochondrial (human convention).
            For other species, you may need to adjust these values based on their
            chromosome structure. The mapper will use species-specific chromosome
            notation but these numbers define the numeric representation.
        log : Log, optional
            Logging object for warnings and messages.
        verbose : bool, default=True
            Whether to print warnings and messages.
        """
        self.species = species.lower() if species else "homo sapiens"
        self.build = build
        self.xymt_num = xymt_num
        self.verbose = verbose
        
        if log is None:
            self.log = Log()
        else:
            self.log = log
        
        # Initialize Chromosomes object for species-specific mappings
        self.chromosomes_obj = Chromosomes(species=self.species)
        
        # Get chromosome mappings
        x, y, mt = self.chromosomes_obj.get_chromosome_mappings(xymt_num=self.xymt_num)
        self.x_mapping = x
        self.y_mapping = y
        self.mt_mapping = mt
        
        # Get species-specific chromosome information
        self.autosomes = self.chromosomes_obj.autosomes
        self.sex_chromosomes = self.chromosomes_obj.sex_chromosomes
        self.mitochondrial = self.chromosomes_obj.mitochondrial
        self.all_chromosomes = self.chromosomes_obj.chromosomes
        
        # Calculate maximum chromosome number for this species
        self._max_chr = self._calculate_max_chromosome()
        
        # Build middle layer mappings (number ↔ string)
        self._build_middle_layer_mappings()
        
        # NCBI RefSeq ID mappings (if build is provided and valid)
        self._nc_mappings = {}
        self._nc_mappings_inverse = {}
        self._num_to_nc = None
        self._nc_to_num = None
        if self.build is not None and self.build not in ("99", "Unknown"):
            self._build_nc_mappings()
        
        # Sumstats layer mappings (will be built when sumstats format is detected)
        self._sumstats_format: Optional[str] = None
        self._sumstats_to_num: Dict[Union[str, int], int] = {}
        self._num_to_sumstats: Dict[int, Union[str, int]] = {}
        self._sumstats_to_num_lower: Dict[str, int] = {}  # Case-insensitive lookup cache
        
        # Reference layer mappings (will be built when reference format is detected)
        self._reference_format: Optional[str] = None
        self._reference_prefix: str = ""
        self._reference_to_num: Dict[Union[str, int], int] = {}
        self._num_to_reference: Dict[int, Union[str, int]] = {}
        self._reference_to_num_lower: Dict[str, int] = {}  # Case-insensitive lookup cache
        
        # Cache for reference file format detection
        self._reference_cache: Dict[str, Tuple[str, str]] = {}
        # Default reference file (set by from_reference_file)
        self._default_reference_file: Optional[str] = None
    
    def _calculate_max_chromosome(self) -> int:
        """Calculate maximum chromosome number for this species."""
        # Get all numeric chromosome values
        max_num = 0
        
        # Check autosomes
        for autosome in self.autosomes:
            try:
                num = int(autosome)
                max_num = max(max_num, num)
            except ValueError:
                pass
        
        # Check sex chromosomes and mitochondrial (from xymt_num)
        if self.xymt_num:
            max_num = max(max_num, max(self.xymt_num))
        
        # Default to at least the highest xymt_num if no autosomes
        if max_num == 0 and self.xymt_num:
            max_num = max(self.xymt_num)
        
        # Ensure minimum of 25 for backward compatibility, but use species-specific if larger
        return max(max_num, 25) if max_num > 0 else 25
    
    def _build_middle_layer_mappings(self):
        """Build mappings for the numeric middle layer using species-specific chromosomes."""
        # String to numeric mapping - use species-specific max_chr
        self._str_to_num = self.chromosomes_obj.get_chr_to_number_dict(
            out_chr=False, xymt_num=self.xymt_num, max_chr=self._max_chr
        )
        # Also support "M" as alias for mitochondrial (species-specific)
        if self.mitochondrial and self.mitochondrial in self._str_to_num:
            # Support common aliases
            mt_aliases = ["M", "m", "mt", "Mt", "MT"]
            mt_num = self._str_to_num[self.mitochondrial]
            for alias in mt_aliases:
                if alias != self.mitochondrial:
                    self._str_to_num[alias] = mt_num
        
        # Numeric to string mapping - use species-specific max_chr
        self._num_to_str = self.chromosomes_obj.get_number_to_chr_dict(
            in_chr=False, xymt_num=self.xymt_num, prefix="", max_chr=self._max_chr
        )
        
        # Fill in any missing numeric keys for the range [1, max_chr]
        for num in range(1, self._max_chr + 1):
            if num not in self._num_to_str:
                if num in self.xymt_num:
                    # Map sex chromosomes and mitochondrial
                    idx = self.xymt_num.index(num)
                    if idx == 0 and len(self.sex_chromosomes) > 0:
                        self._num_to_str[num] = self.x_mapping[0]
                    elif idx == 1 and len(self.sex_chromosomes) > 1:
                        self._num_to_str[num] = self.y_mapping[0]
                    elif idx == 2 and self.mitochondrial:
                        self._num_to_str[num] = self.mt_mapping[0]
                    else:
                        self._num_to_str[num] = str(num)
                else:
                    # Regular autosome
                    self._num_to_str[num] = str(num)
    
    def _build_nc_mappings(self):
        """Build NCBI RefSeq ID mappings."""
        if self.build is None:
            return
        
        # Get string to NC mapping
        chr_to_nc = get_chr_to_NC(
            build=self.build,
            inverse=False,
            species=self.species,
            log=self.log,
            verbose=self.verbose
        )
        
        if not chr_to_nc:
            return
        
        # Build NC to string mapping
        self._nc_mappings = chr_to_nc.copy()
        self._nc_mappings_inverse = {v: k for k, v in chr_to_nc.items()}
        
        # Build numeric to NC mapping
        self._num_to_nc = {}
        self._nc_to_num = {}
        for chrom_str, nc_id in chr_to_nc.items():
            if chrom_str in self._str_to_num:
                chrom_num = self._str_to_num[chrom_str]
                self._num_to_nc[chrom_num] = nc_id
                self._nc_to_num[nc_id] = chrom_num
    
    # ============================================================================
    # Format Detection
    # ============================================================================
    
    def detect_format(self, data: Union[pd.Series, list, Any], sample_size: int = 100) -> str:
        """
        Detect chromosome format from data (species-aware).
        
        Parameters
        ----------
        data : pd.Series, list, or single value
            Chromosome data to analyze.
        sample_size : int, default=100
            Number of samples to check for format detection.
        
        Returns
        -------
        str
            Detected format: "nc", "chr", "numeric", or "string".
        """
        # Handle single value
        if not isinstance(data, (pd.Series, list)):
            if isinstance(data, (int, np.integer)):
                return "numeric"
            data_str = str(data).strip()
            if self._NCBI_PATTERN.match(data_str):
                return "nc"
            if self._CHR_PREFIX_PATTERN.match(data_str):
                return "chr"
            # Check if it matches species-specific sex chromosomes
            data_upper = data_str.upper()
            if data_upper in [c.upper() for c in self.sex_chromosomes]:
                return "string"
            if self.mitochondrial and data_upper == self.mitochondrial.upper():
                return "string"
            return "string"
        
        # Handle Series or list
        if isinstance(data, pd.Series):
            series = data.dropna()
        else:
            series = pd.Series(data).dropna()
        
        if len(series) == 0:
            return "string"  # Default
        
        # Sample for efficiency
        sample = series.head(min(sample_size, len(series)))
        
        # Check for NCBI RefSeq IDs
        sample_str = sample.astype(str).str.strip()
        # Convert compiled pattern to string for pandas str.match
        nc_pattern_str = self._NCBI_PATTERN.pattern
        if sample_str.str.match(nc_pattern_str, case=False, na=False).any():
            return "nc"
        
        # Check for chr prefix
        chr_pattern_str = self._CHR_PREFIX_PATTERN.pattern
        if sample_str.str.match(chr_pattern_str, case=False, na=False).any():
            return "chr"
        
        # Check if values match species-specific chromosomes (case-insensitive)
        # This helps identify string format, but we continue to numeric check
        series_upper = sample_str.str.upper()
        for chrom in self.all_chromosomes:
            if series_upper.str.contains(f"^{re.escape(chrom.upper())}$", case=False, na=False).any():
                # Found species-specific chromosome, likely string format
                break
        
        # Check if numeric
        if pd.api.types.is_integer_dtype(sample):
            return "numeric"
        
        # Try to convert to numeric
        try:
            numeric_sample = pd.to_numeric(sample_str, errors='coerce')
            if numeric_sample.notna().all():
                return "numeric"
        except (ValueError, TypeError):
            pass
        
        # Default to string
        return "string"
    
    def detect_sumstats_format(self, data: Union[pd.Series, list, Any], sample_size: int = 100) -> str:
        """
        Detect sumstats chromosome format and build sumstats layer mappings.
        
        Parameters
        ----------
        data : pd.Series, list, or single value
            Sumstats chromosome data to analyze.
        sample_size : int, default=100
            Number of samples to check for format detection.
        
        Returns
        -------
        str
            Detected format: "nc", "chr", "numeric", or "string".
        """
        format_type = self.detect_format(data, sample_size)
        self._build_sumstats_layer(format_type)
        return format_type
    
    def detect_reference_format(self, reference_file: str) -> Tuple[str, str]:
        """
        Detect reference file chromosome format and build reference layer mappings.
        
        Parameters
        ----------
        reference_file : str
            Path to reference file (VCF, GTF, FASTA, chain file).
        
        Returns
        -------
        tuple[str, str]
            Tuple of (format_type, prefix) where format_type is "nc", "chr", "numeric", or "string",
            and prefix is the chr prefix (e.g., "chr", "Chr", "") if format is "chr".
        """
        # Check cache
        if reference_file in self._reference_cache:
            format_type, prefix = self._reference_cache[reference_file]
            self._build_reference_layer(format_type, prefix, reference_file=reference_file)
            return format_type, prefix
        
        # Detect file type and sample chromosomes
        format_type, prefix = self._detect_reference_format_from_file(reference_file)
        
        # Cache result
        self._reference_cache[reference_file] = (format_type, prefix)
        
        # Build reference layer (pass reference_file for NCBI format to extract contigs)
        self._build_reference_layer(format_type, prefix, reference_file=reference_file)
        
        return format_type, prefix
    
    def _detect_reference_format_from_file(self, reference_file: str) -> Tuple[str, str]:
        """Detect chromosome format from reference file."""
        if not os.path.exists(reference_file):
            if self.verbose:
                self.log.warning(f"Reference file not found: {reference_file}. Using default format.")
            return "string", ""
        
        # Lazy imports to avoid circular dependency
        from gwaslab.io.io_vcf import VCF_BCF_SUFFIXES
        from gwaslab.io.io_gtf import GTF_GFF_SUFFIXES
        from gwaslab.io.io_fasta import FASTA_SUFFIXES
        from gwaslab.io.io_chain import CHAIN_SUFFIXES
        
        # Sample chromosomes from file based on file type
        # Handle double extensions like .vcf.gz, .gtf.gz, etc.
        file_ext_lower = reference_file.lower()
        
        # Optimized file type detection using tuple checks
        if file_ext_lower.endswith(VCF_BCF_SUFFIXES):
            return self._detect_vcf_format(reference_file)
        elif file_ext_lower.endswith(GTF_GFF_SUFFIXES):
            return self._detect_gtf_format(reference_file)
        elif file_ext_lower.endswith(FASTA_SUFFIXES):
            return self._detect_fasta_format(reference_file)
        elif file_ext_lower.endswith(CHAIN_SUFFIXES):
            return self._detect_chain_format(reference_file)
        else:
            # Default: try to detect from first few lines
            return self._detect_generic_format(reference_file)
    
    def _extract_contigs_from_file(self, reference_file: str) -> list:
        """Extract contig/chromosome names from reference file."""
        # Lazy imports to avoid circular dependency
        from gwaslab.io.io_vcf import VCF_BCF_SUFFIXES
        from gwaslab.io.io_gtf import GTF_GFF_SUFFIXES
        from gwaslab.io.io_fasta import FASTA_SUFFIXES
        
        contigs = []
        file_ext_lower = reference_file.lower()
        
        if file_ext_lower.endswith(VCF_BCF_SUFFIXES):
            try:
                import pysam
                with pysam.VariantFile(reference_file) as vcf:
                    contigs = list(vcf.header.contigs)
            except Exception as e:
                if self.verbose:
                    self.log.warning(f"Could not extract contigs from VCF: {e}")
                pass
        elif file_ext_lower.endswith(GTF_GFF_SUFFIXES):
            try:
                import polars as pl
                df = pl.read_csv(reference_file, separator='\t', comment_char='#', has_header=False, n_rows=1000)
                if df.shape[1] > 0:
                    contigs = df.select(df.columns[0]).unique().to_series().to_list()
            except (ImportError, Exception):
                pass
        elif file_ext_lower.endswith(FASTA_SUFFIXES):
            try:
                import pysam
                with pysam.FastaFile(reference_file) as fasta:
                    contigs = list(fasta.references)
            except (ImportError, Exception):
                pass
        
        return contigs
    
    def _detect_vcf_format(self, vcf_path: str) -> Tuple[str, str]:
        """Detect chromosome format from VCF file."""
        try:
            import pysam
            with pysam.VariantFile(vcf_path) as vcf:
                # First try to detect from header contigs (more reliable)
                contigs = list(vcf.header.contigs)
                if contigs:
                    # Sample up to 100 contigs
                    sample_contigs = contigs[:min(100, len(contigs))]
                    format_type = self.detect_format(sample_contigs)
                    prefix = self._extract_prefix(sample_contigs[0]) if format_type == "chr" else ""
                    if format_type in ["nc", "chr", "string"]:
                        return format_type, prefix
                
                # Fallback: read from records
                chromosomes = []
                count = 0
                for record in vcf:
                    chrom = record.chrom
                    chromosomes.append(str(chrom))
                    count += 1
                    if count >= 100:
                        break
                
                if chromosomes:
                    format_type = self.detect_format(chromosomes)
                    prefix = self._extract_prefix(chromosomes[0]) if format_type == "chr" else ""
                    return format_type, prefix
        except Exception as e:
            if self.verbose:
                self.log.warning(f"Could not read VCF file with pysam: {e}. Trying alternative method.")
        
        # Fallback: read as text
        try:
            import gzip
            open_func = gzip.open if vcf_path.endswith('.gz') else open
            mode = 'rt' if vcf_path.endswith('.gz') else 'r'
            with open_func(vcf_path, mode) as f:
                chromosomes = []
                for line in f:
                    if line.startswith('#CHROM'):
                        # Header line, skip
                        continue
                    if line.startswith('#'):
                        continue
                    parts = line.split('\t')
                    if len(parts) > 0:
                        chromosomes.append(parts[0])
                        if len(chromosomes) >= 100:
                            break
                
                if chromosomes:
                    format_type = self.detect_format(chromosomes)
                    prefix = self._extract_prefix(chromosomes[0]) if format_type == "chr" else ""
                    return format_type, prefix
        except Exception as e:
            if self.verbose:
                self.log.warning(f"Could not read VCF file: {e}.")
        
        return "string", ""
    
    def _detect_gtf_format(self, gtf_path: str) -> Tuple[str, str]:
        """Detect chromosome format from GTF file."""
        try:
            import gzip
            open_func = gzip.open if gtf_path.endswith('.gz') else open
            mode = 'rt' if gtf_path.endswith('.gz') else 'r'
            with open_func(gtf_path, mode) as f:
                chromosomes = []
                for line in f:
                    if line.startswith('#'):
                        continue
                    parts = line.split('\t')
                    if len(parts) > 0:
                        chromosomes.append(parts[0])
                        if len(chromosomes) >= 100:
                            break
                
                if chromosomes:
                    format_type = self.detect_format(chromosomes)
                    prefix = self._extract_prefix(chromosomes[0]) if format_type == "chr" else ""
                    return format_type, prefix
        except Exception as e:
            if self.verbose:
                self.log.warning(f"Could not read GTF file: {e}.")
        
        return "string", ""
    
    def _detect_fasta_format(self, fasta_path: str) -> Tuple[str, str]:
        """Detect chromosome format from FASTA file."""
        try:
            import gzip
            open_func = gzip.open if fasta_path.endswith('.gz') or fasta_path.endswith('.bgz') else open
            mode = 'rt' if (fasta_path.endswith('.gz') or fasta_path.endswith('.bgz')) else 'r'
            with open_func(fasta_path, mode) as f:
                chromosomes = []
                for line in f:
                    if line.startswith('>'):
                        chrom = line[1:].strip().split()[0]  # Get first word after >
                        chromosomes.append(chrom)
                        if len(chromosomes) >= 100:
                            break
                
                if chromosomes:
                    format_type = self.detect_format(chromosomes)
                    prefix = self._extract_prefix(chromosomes[0]) if format_type == "chr" else ""
                    return format_type, prefix
        except Exception as e:
            if self.verbose:
                self.log.warning(f"Could not read FASTA file: {e}.")
        
        return "string", ""
    
    def _detect_chain_format(self, chain_path: str) -> Tuple[str, str]:
        """Detect chromosome format from chain file."""
        try:
            import gzip
            open_func = gzip.open if chain_path.endswith('.gz') else open
            mode = 'rt' if chain_path.endswith('.gz') else 'r'
            with open_func(chain_path, mode) as f:
                chromosomes = []
                for line in f:
                    if line.startswith('chain'):
                        parts = line.split()
                        if len(parts) >= 3:
                            chromosomes.append(parts[2])  # Target chromosome
                            if len(chromosomes) >= 100:
                                break
                
                if chromosomes:
                    format_type = self.detect_format(chromosomes)
                    prefix = self._extract_prefix(chromosomes[0]) if format_type == "chr" else ""
                    return format_type, prefix
        except Exception as e:
            if self.verbose:
                self.log.warning(f"Could not read chain file: {e}.")
        
        return "string", ""
    
    def _detect_generic_format(self, file_path: str) -> Tuple[str, str]:
        """Generic format detection from file."""
        try:
            import gzip
            open_func = gzip.open if file_path.endswith('.gz') else open
            mode = 'rt' if file_path.endswith('.gz') else 'r'
            with open_func(file_path, mode) as f:
                chromosomes = []
                for i, line in enumerate(f):
                    if i >= 100:
                        break
                    # Try to extract chromosome from first column
                    parts = line.strip().split()
                    if parts:
                        chromosomes.append(parts[0])
                
                if chromosomes:
                    format_type = self.detect_format(chromosomes)
                    prefix = self._extract_prefix(chromosomes[0]) if format_type == "chr" else ""
                    return format_type, prefix
        except Exception as e:
            if self.verbose:
                self.log.warning(f"Could not read file: {e}.")
        
        return "string", ""
    
    def _extract_prefix(self, chrom: str) -> str:
        """Extract chr prefix from chromosome string."""
        match = self._CHR_PREFIX_PATTERN.match(chrom)
        if match:
            return match.group(1)
        return ""
    
    # ============================================================================
    # Layer Building
    # ============================================================================
    
    def _build_sumstats_layer(self, format_type: str):
        """Build sumstats layer mappings based on format type using species-specific chromosomes."""
        self._sumstats_format = format_type
        self._sumstats_to_num = {}
        self._num_to_sumstats = {}
        self._sumstats_to_num_lower = {}  # Clear case-insensitive cache
        
        if format_type == "numeric":
            # Direct numeric mapping - use species-specific max_chr
            # Use dict comprehension for better performance
            self._sumstats_to_num = {num: num for num in range(1, self._max_chr + 1)}
            self._num_to_sumstats = {num: num for num in range(1, self._max_chr + 1)}
        elif format_type == "string":
            # String to number mapping - use dict comprehension
            self._sumstats_to_num = dict(self._str_to_num)
            self._num_to_sumstats = {num: chrom_str for chrom_str, num in self._str_to_num.items()}
        elif format_type == "chr":
            # Chr-prefixed to number mapping
            # Pre-compute prefixes once
            prefixes = ["chr", "Chr", "CHR"]
            for chrom_str, num in self._str_to_num.items():
                # Build chr-prefixed versions
                for prefix in prefixes:
                    chr_key = prefix + chrom_str
                    self._sumstats_to_num[chr_key] = num
                    self._sumstats_to_num[chr_key.lower()] = num
                # Store reverse mapping (use first prefix for consistency)
                self._num_to_sumstats[num] = "chr" + chrom_str
        elif format_type == "nc":
            # NCBI RefSeq ID to number mapping
            if self._nc_to_num:
                self._sumstats_to_num = dict(self._nc_to_num)
                self._num_to_sumstats = {num: nc_id for nc_id, num in self._nc_to_num.items()}
    
    def _build_reference_layer(self, format_type: str, prefix: str = "", reference_file: Optional[str] = None):
        """Build reference layer mappings based on format type and prefix using species-specific chromosomes."""
        self._reference_format = format_type
        self._reference_prefix = prefix
        self._reference_to_num = {}
        self._num_to_reference = {}
        self._reference_to_num_lower = {}  # Clear case-insensitive cache
        
        # Dispatch to format-specific builders
        builders = {
            "numeric": self._build_reference_numeric,
            "string": self._build_reference_string,
            "chr": self._build_reference_chr,
            "nc": self._build_reference_nc,
        }
        
        builder = builders.get(format_type)
        if builder:
            builder(prefix, reference_file)
        else:
            if self.verbose:
                self.log.warning(f"Unknown format type: {format_type}. Using default string format.")
            self._build_reference_string(prefix, reference_file)
    
    def _build_reference_numeric(self, prefix: str, reference_file: Optional[str]):
        """Build numeric format reference layer."""
        # Direct numeric mapping - use species-specific max_chr
        # Use dict comprehension for better performance
        self._reference_to_num = {num: num for num in range(1, self._max_chr + 1)}
        self._num_to_reference = {num: num for num in range(1, self._max_chr + 1)}
    
    def _build_reference_string(self, prefix: str, reference_file: Optional[str]):
        """Build string format reference layer."""
        # String to number mapping - use dict comprehension
        self._reference_to_num = dict(self._str_to_num)
        self._num_to_reference = {num: chrom_str for chrom_str, num in self._str_to_num.items()}
    
    def _build_reference_chr(self, prefix: str, reference_file: Optional[str]):
        """Build chr-prefixed format reference layer."""
        # Use provided prefix or default to "chr"
        actual_prefix = prefix if prefix else "chr"
        for chrom_str, num in self._str_to_num.items():
            chr_key = actual_prefix + chrom_str
            self._reference_to_num[chr_key] = num
            self._reference_to_num[chr_key.lower()] = num
            # Store reverse mapping
            self._num_to_reference[num] = actual_prefix + chrom_str
    
    def _build_reference_nc(self, prefix: str, reference_file: Optional[str]):
        """Build NCBI RefSeq ID format reference layer."""
        if self._nc_to_num:
            # Use pre-built mappings if available (requires build)
            self._reference_to_num = dict(self._nc_to_num)
            self._num_to_reference = {num: nc_id for nc_id, num in self._nc_to_num.items()}
        else:
            # Build mappings from actual contigs if build not provided
            self._build_reference_nc_from_file(reference_file)
    
    def _build_reference_nc_from_file(self, reference_file: Optional[str]):
        """Build NCBI mappings from reference file contigs.
        
        If build is not set, this method will attempt to auto-detect the build
        by matching contigs against known NCBI RefSeq IDs for different builds.
        """
        if not reference_file or not os.path.exists(reference_file):
            if self.verbose:
                self.log.warning(
                    "NCBI format detected but could not build mappings. "
                    "Specify build parameter or ensure reference file is accessible."
                )
            return
        
        try:
            contigs = self._extract_contigs_from_file(reference_file)
            if not contigs:
                if self.verbose:
                    self.log.warning(
                        "NCBI format detected but no contigs found in reference file. "
                        "Specify build parameter or ensure reference file is accessible."
                    )
                return
            
            # If build is not set, or is "99" (unknown), or NC mappings are not available,
            # try to auto-detect it from the contigs
            build_is_invalid = (self.build is None or 
                               self.build == "99" or 
                               self.build == "Unknown" or
                               not self._nc_to_num)
            
            if build_is_invalid:
                detected_build = self._detect_build_from_nc_contigs(contigs)
                if detected_build:
                    if self.verbose:
                        self.log.write(f" -Auto-detected build {detected_build} from NC notation in reference file...", verbose=self.verbose)
                    self.build = detected_build
                    self._build_nc_mappings()
                    # Now use the pre-built mappings
                    self._reference_to_num = dict(self._nc_to_num)
                    self._num_to_reference = {num: nc_id for nc_id, num in self._nc_to_num.items()}
                    return
                else:
                    if self.verbose:
                        self.log.warning(
                            "NCBI format detected but could not auto-detect build. "
                            "Specify build parameter for proper NC mappings."
                        )
                    return
            
            # If build is set and valid, use the pre-built mappings
            if self._nc_to_num:
                self._reference_to_num = dict(self._nc_to_num)
                self._num_to_reference = {num: nc_id for nc_id, num in self._nc_to_num.items()}
            else:
                # Fallback: parse NCBI IDs directly (may not work for all species)
                for nc_id in contigs:
                    chrom_num = self._parse_ncbi_id_to_number(nc_id)
                    if chrom_num is not None:
                        self._reference_to_num[nc_id] = chrom_num
                        # Only set reverse mapping if not already set (prefer first seen)
                        if chrom_num not in self._num_to_reference:
                            self._num_to_reference[chrom_num] = nc_id
        except Exception as e:
            if self.verbose:
                self.log.warning(f"Could not extract contigs from {reference_file}: {e}")
    
    def _detect_build_from_nc_contigs(self, contigs: list) -> Optional[str]:
        """Auto-detect genome build from NCBI RefSeq contig IDs.
        
        Parameters
        ----------
        contigs : list
            List of contig/chromosome names from the reference file.
        
        Returns
        -------
        str or None
            Detected build (e.g., "19", "38") if found, None otherwise.
        """
        from gwaslab.bd.bd_common_data import get_number_to_NC
        
        # Check which build matches
        for build_candidate in ["19", "38"]:
            nc_mapping = get_number_to_NC(build=build_candidate, species=self.species, log=self.log, verbose=False)
            if nc_mapping:
                nc_values = set(nc_mapping.values())
                # Check if any contig matches this build's NC IDs
                if any(contig in nc_values for contig in contigs):
                    return build_candidate
        
        return None
    
    def _parse_ncbi_id_to_number(self, nc_id: str) -> Optional[int]:
        """Parse NCBI RefSeq ID to chromosome number using species-specific mapping.
        
        This function uses the species-specific NCBI mapping (from _nc_to_num)
        to correctly convert NCBI RefSeq IDs to chromosome numbers for all species.
        Direct parsing is not used because the numeric part of NCBI IDs does not
        correspond to chromosome numbers for non-human species. For example:
        - Human: NC_000001.10 -> 1 (numeric part happens to match)
        - Mouse: NC_000067.7 -> 1 (numeric part is 67, but chromosome is 1)
        
        Parameters
        ----------
        nc_id : str
            NCBI RefSeq ID (e.g., "NC_000001.10", "NC_000067.7")
        
        Returns
        -------
        int or None
            Chromosome number if valid mapping exists, None otherwise
        """
        # Validate format with regex
        if not self._NCBI_PATTERN.match(nc_id):
            return None
        
        # Use species-specific mapping (required for correct conversion across all species)
        if self._nc_to_num and nc_id in self._nc_to_num:
            chrom_num = self._nc_to_num[nc_id]
            if 1 <= chrom_num <= self._max_chr:
                return chrom_num
        
        # No mapping available - return None
        # This ensures we don't incorrectly parse IDs for non-human species
        return None
    
    # ============================================================================
    # Core Mapping Methods: Middle Layer
    # ============================================================================
    
    def _to_number(self, chromosome: Union[str, int], layer: str = "sumstats") -> Union[int, float]:
        """
        Convert chromosome to number (middle layer).
        
        Parameters
        ----------
        chromosome : str or int
            Chromosome identifier.
        layer : str, default="sumstats"
            Which layer to use: "sumstats" or "reference".
        
        Returns
        -------
        int or float
            Chromosome number, or pd.NA if conversion fails.
        """
        if isinstance(chromosome, (int, np.integer)):
            return int(chromosome)
        
        chrom_str = str(chromosome).strip()
        
        # Use appropriate layer mapping
        if layer == "sumstats":
            mapping = self._sumstats_to_num
        else:
            mapping = self._reference_to_num
        
        # Try direct lookup
        if chrom_str in mapping:
            return mapping[chrom_str]
        
        # Try case-insensitive lookup using cache
        chrom_lower = chrom_str.lower()
        if layer == "sumstats":
            lower_cache = self._sumstats_to_num_lower
        else:
            lower_cache = self._reference_to_num_lower
        
        # Build cache if needed (lazy initialization)
        # Only build if cache is empty but mapping has items
        if len(lower_cache) == 0 and len(mapping) > 0:
            for key, value in mapping.items():
                if isinstance(key, str):
                    lower_cache[key.lower()] = value
        
        # Use cache for fast lookup
        if chrom_lower in lower_cache:
            return lower_cache[chrom_lower]
        
        # Try to parse as number
        try:
            num = int(chrom_str)
            if 1 <= num <= self._max_chr:
                return num
        except (ValueError, TypeError):
            pass
        
        # Try NCBI RefSeq ID
        if self._NCBI_PATTERN.match(chrom_str):
            if chrom_str in self._nc_to_num:
                return self._nc_to_num[chrom_str]
        
        # Default: try string to number mapping
        if chrom_str in self._str_to_num:
            return self._str_to_num[chrom_str]
        
        # If all else fails, return pd.NA instead of raising error
        return pd.NA
    
    def _from_number(self, number: int, layer: str = "sumstats") -> Union[str, int]:
        """
        Convert number (middle layer) to chromosome format.
        
        Parameters
        ----------
        number : int
            Chromosome number.
        layer : str, default="sumstats"
            Which layer to use: "sumstats" or "reference".
        
        Returns
        -------
        str or int
            Chromosome identifier in the specified layer format.
        """
        if not (1 <= number <= self._max_chr):
            raise ValueError(
                f"Invalid chromosome number for {self.species}: {number}. "
                f"Valid range: 1-{self._max_chr}"
            )
        
        # Use appropriate layer mapping
        if layer == "sumstats":
            mapping = self._num_to_sumstats
        else:
            mapping = self._num_to_reference
        
        if number in mapping:
            return mapping[number]
        
        # Fallback to middle layer
        if number in self._num_to_str:
            return self._num_to_str[number]
        
        # Ultimate fallback
        return number
    
    # ============================================================================
    # Public API: Main Mapping Methods
    # ============================================================================
    
    def sumstats_to_number(self, chromosome: Union[str, int]) -> Union[int, float]:
        """
        Convert sumstats chromosome to number (middle layer).
        
        Parameters
        ----------
        chromosome : str or int
            Sumstats chromosome identifier.
        
        Returns
        -------
        int or float
            Chromosome number, or pd.NA if conversion fails.
        
        Notes
        -----
        If `detect_sumstats_format()` has not been called, this method will:
        - Default to numeric format (assumes sumstats uses numeric chromosomes)
        - Work for integers (1, 2, 23) - returns as-is
        - Work for numeric strings ("1", "2", "23") - parses as int if in valid range
        - Work for string chromosomes ("X", "Y", "MT") - uses fallback string mapping
        - For format-specific strings like "chr1", will attempt conversion using fallback logic
        
        It is recommended to call `detect_sumstats_format()` first to ensure
        all format-specific conversions work correctly.
        """
        # If sumstats layer not built, default to numeric format
        if self._sumstats_format is None:
            self._build_sumstats_layer("numeric")
        return self._to_number(chromosome, layer="sumstats")
    
    def number_to_sumstats(self, number: int) -> Union[str, int]:
        """
        Convert number (middle layer) to sumstats chromosome format.
        
        Parameters
        ----------
        number : int
            Chromosome number.
        
        Returns
        -------
        str or int
            Sumstats chromosome identifier.
        
        Notes
        -----
        If `detect_sumstats_format()` has not been called, this method will:
        - Default to numeric format (returns the number itself)
        - This ensures consistent behavior when sumstats format is unknown
        
        It is recommended to call `detect_sumstats_format()` first to ensure
        the output matches the actual sumstats format.
        """
        # If sumstats layer not built, default to numeric format
        if self._sumstats_format is None:
            self._build_sumstats_layer("numeric")
        return self._from_number(number, layer="sumstats")
    
    def reference_to_number(self, chromosome: Union[str, int], reference_file: Optional[str] = None) -> Union[int, float]:
        """
        Convert reference chromosome to number (middle layer).
        
        Parameters
        ----------
        chromosome : str or int
            Reference chromosome identifier.
        reference_file : str, optional
            Path to reference file. If provided and reference layer not yet built,
            will detect format from file.
        
        Returns
        -------
        int or float
            Chromosome number, or pd.NA if conversion fails.
        """
        if reference_file is not None:
            self.detect_reference_format(reference_file)
        return self._to_number(chromosome, layer="reference")
    
    def number_to_reference(self, number: int, reference_file: Optional[str] = None) -> Union[str, int]:
        """
        Convert number (middle layer) to reference chromosome format.
        
        Parameters
        ----------
        number : int
            Chromosome number.
        reference_file : str, optional
            Path to reference file. If provided and reference layer not yet built,
            will detect format from file.
        
        Returns
        -------
        str or int
            Reference chromosome identifier.
        """
        if reference_file is not None:
            self.detect_reference_format(reference_file)
        return self._from_number(number, layer="reference")
    
    def sumstats_to_reference(self, 
                             chromosome: Union[str, int], 
                             reference_file: Optional[str] = None,
                             as_string: bool = False) -> Union[str, int]:
        """
        Convert sumstats chromosome to reference chromosome format.
        
        This is the main method for matching sumstats chr with reference chr.
        
        Parameters
        ----------
        chromosome : str or int
            Sumstats chromosome identifier.
        reference_file : str, optional
            Path to reference file. If provided and reference layer not yet built,
            will detect format from file.
        as_string : bool, default=False
            If True, always return string (useful for VCF compatibility with pysam).
            If False, return the original type (str or int) as determined by reference format.
        
        Returns
        -------
        str or int
            Reference chromosome identifier.
            Returns string representation of pd.NA if conversion fails.
        
        Examples
        --------
        >>> mapper = ChromosomeMapper()
        >>> mapper.detect_sumstats_format(pd.Series([1, 2, 23]))
        'numeric'
        >>> mapper.sumstats_to_reference(1, reference_file="reference.vcf.gz")
        'chr1'
        >>> mapper.sumstats_to_reference(1, reference_file="reference.vcf.gz", as_string=True)
        'chr1'
        """
        # Convert sumstats → number → reference
        # sumstats_to_number will auto-default to numeric if layer not built
        number = self.sumstats_to_number(chromosome)
        if reference_file is not None:
            self.detect_reference_format(reference_file)
        result = self.number_to_reference(number, reference_file)
        
        # Convert to string if requested (for VCF compatibility with pysam)
        if as_string:
            if pd.isna(result):
                return str(result)  # Returns "NA" for pd.NA
            return str(result)
        
        # Return original type
        return result
    
    def reference_to_sumstats(self, 
                             chromosome: Union[str, int], 
                             reference_file: Optional[str] = None) -> Union[str, int]:
        """
        Convert reference chromosome to sumstats chromosome format.
        
        This is the main method for matching reference chr with sumstats chr.
        
        Parameters
        ----------
        chromosome : str or int
            Reference chromosome identifier.
        reference_file : str, optional
            Path to reference file. If provided and reference layer not yet built,
            will detect format from file.
        
        Returns
        -------
        str or int
            Sumstats chromosome identifier.
        
        Examples
        --------
        >>> mapper = ChromosomeMapper()
        >>> mapper.detect_sumstats_format(pd.Series([1, 2, 23]))
        'numeric'
        >>> mapper.reference_to_sumstats("chr1", reference_file="reference.vcf.gz")
        1
        """
        # Convert reference → number → sumstats
        if reference_file is not None:
            self.detect_reference_format(reference_file)
        number = self.reference_to_number(chromosome, reference_file)
        return self.number_to_sumstats(number)
    
    # ============================================================================
    # Convenience Methods: Series and DataFrame Operations
    # ============================================================================
    
    def sumstats_to_reference_series(self, 
                                    series: pd.Series, 
                                    reference_file: Optional[str] = None) -> pd.Series:
        """
        Convert sumstats chromosome Series to reference format.
        
        Parameters
        ----------
        series : pd.Series
            Sumstats chromosome Series.
        reference_file : str, optional
            Path to reference file.
        
        Returns
        -------
        pd.Series
            Reference chromosome Series.
        """
        if reference_file is not None:
            self.detect_reference_format(reference_file)
        
        # Ensure sumstats layer is built
        if self._sumstats_format is None:
            self._build_sumstats_layer("numeric")
        
        # Ensure reference layer is built
        if self._reference_format is None:
            if reference_file is not None:
                self.detect_reference_format(reference_file)
            else:
                # Default to string format if no reference file
                self._build_reference_layer("string")
        
        # Vectorized conversion: sumstats -> number -> reference
        # Convert to number first (vectorized using map)
        num_series = series.map(lambda x: self._to_number(x, layer="sumstats"))
        
        # Convert number to reference format (vectorized using map)
        result = num_series.map(lambda x: self._from_number(x, layer="reference"))
        
        # Convert to object dtype to avoid dtype incompatibility warnings
        result = result.astype(object)
        
        return result
    
    def reference_to_sumstats_series(self, 
                                    series: pd.Series, 
                                    reference_file: Optional[str] = None) -> pd.Series:
        """
        Convert reference chromosome Series to sumstats format.
        
        Parameters
        ----------
        series : pd.Series
            Reference chromosome Series.
        reference_file : str, optional
            Path to reference file.
        
        Returns
        -------
        pd.Series
            Sumstats chromosome Series.
        """
        if reference_file is not None:
            self.detect_reference_format(reference_file)
        
        # Convert to object dtype to avoid dtype incompatibility warnings
        result = series.copy().astype(object)
        for idx, chrom in series.items():
            try:
                result.loc[idx] = self.reference_to_sumstats(chrom, reference_file)
            except (ValueError, KeyError):
                result.loc[idx] = chrom  # Keep original if conversion fails
        
        return result
    
    # ============================================================================
    # Class Methods: Factory Methods
    # ============================================================================
    
    @classmethod
    def from_reference_file(cls,
                           file_path: str,
                           species: str = "homo sapiens",
                           build: Optional[str] = None,
                           log: Optional[Log] = None,
                           verbose: bool = True) -> "ChromosomeMapper":
        """
        Create a ChromosomeMapper and detect reference format from file.
        
        Parameters
        ----------
        file_path : str
            Path to reference file.
        species : str, default="homo sapiens"
            Species name.
        build : str, optional
            Genome build version.
        log : Log, optional
            Logging object.
        verbose : bool, default=True
            Whether to print messages.
        
        Returns
        -------
        ChromosomeMapper
            Mapper instance with reference format detected.
        """
        mapper = cls(species=species, build=build, log=log, verbose=verbose)
        mapper._default_reference_file = file_path
        mapper.detect_reference_format(file_path)
        return mapper
    
    # ============================================================================
    # Backward Compatibility Methods
    # ============================================================================
    
    def map(self, chromosome: Union[str, int], reference_file: Optional[str] = None) -> Union[str, int]:
        """
        Map chromosome (backward compatibility method).
        
        If reference_file is provided, maps sumstats format to reference format.
        Otherwise, returns in sumstats format.
        """
        if reference_file is not None:
            return self.sumstats_to_reference(chromosome, reference_file)
        else:
            return self.number_to_sumstats(self.sumstats_to_number(chromosome))
    
    def map_series(self, series: pd.Series, reference_file: Optional[str] = None) -> pd.Series:
        """
        Map chromosome Series (backward compatibility method).
        """
        if reference_file is not None:
            return self.sumstats_to_reference_series(series, reference_file)
        else:
            return series.apply(lambda x: self.number_to_sumstats(self.sumstats_to_number(x)))
    
    def to_numeric(self, chromosome: Union[str, int, pd.Series]) -> Union[int, float, pd.Series]:
        """Convert to numeric format (backward compatibility)."""
        if isinstance(chromosome, pd.Series):
            return chromosome.apply(self.sumstats_to_number)
        return self.sumstats_to_number(chromosome)
    
    def to_string(self, chromosome: Union[str, int, pd.Series]) -> Union[str, pd.Series]:
        """Convert to string format (backward compatibility)."""
        # Convert to number first (will default to numeric if sumstats layer not built)
        if isinstance(chromosome, pd.Series):
            numbers = chromosome.apply(lambda x: self.sumstats_to_number(x))
            # Convert number to string format using _num_to_str mapping
            return numbers.apply(lambda n: self._num_to_str.get(n, str(n)))
        number = self.sumstats_to_number(chromosome)
        # Convert number to string format using _num_to_str mapping
        return self._num_to_str.get(number, str(number))
    
    def to_chr(self, chromosome: Union[str, int, pd.Series], prefix: Optional[str] = None) -> Union[str, pd.Series]:
        """Convert to chr-prefixed format (backward compatibility)."""
        if isinstance(chromosome, pd.Series):
            return chromosome.apply(lambda x: self._to_chr_single(x, prefix))
        return self._to_chr_single(chromosome, prefix)
    
    def _to_chr_single(self, chromosome: Union[str, int], prefix: Optional[str] = None) -> str:
        """Convert single chromosome to chr format."""
        number = self.sumstats_to_number(chromosome)
        chrom_str = self._num_to_str.get(number, str(number))
        actual_prefix = prefix if prefix else "chr"
        return actual_prefix + chrom_str
    
    def to_nc(self, chromosome: Union[str, int, pd.Series]) -> Union[str, pd.Series]:
        """Convert to NCBI RefSeq ID format (backward compatibility)."""
        if isinstance(chromosome, pd.Series):
            return chromosome.apply(lambda x: self._to_nc_single(x))
        return self._to_nc_single(chromosome)
    
    def _to_nc_single(self, chromosome: Union[str, int]) -> str:
        """Convert single chromosome to NC format."""
        number = self.sumstats_to_number(chromosome)
        if not hasattr(self, '_num_to_nc') or self._num_to_nc is None:
            raise ValueError(f"NCBI RefSeq ID mapping not available. Build must be specified.")
        if number in self._num_to_nc:
            return self._num_to_nc[number]
        raise ValueError(f"NCBI RefSeq ID not available for chromosome {chromosome}")
    
    def print(self):
        """Print mapper configuration and state information.
        
        Displays:
        - Species and build information
        - Detected sumstats and reference formats
        - Chromosome mappings
        - Example conversions
        """
        print("=" * 70)
        print("ChromosomeMapper Configuration")
        print("=" * 70)
        
        # Basic information
        print(f"Species: {self.species}")
        print(f"Build: {self.build if self.build else 'Not specified'}")
        print(f"Max chromosome number: {self._max_chr}")
        print(f"X/Y/MT numeric values: {self.xymt_num}")
        print()
        
        # Sumstats format
        if self._sumstats_format:
            print(f"Sumstats format: {self._sumstats_format}")
            if self._sumstats_to_num:
                # Show some example mappings
                examples = list(self._sumstats_to_num.items())[:5]
                print(f"  Example sumstats mappings: {dict(examples)}")
        else:
            print("Sumstats format: Not detected (will default to numeric)")
        print()
        
        # Reference format
        if self._reference_format:
            print(f"Reference format: {self._reference_format}")
            if self._reference_prefix:
                print(f"  Prefix: '{self._reference_prefix}'")
            if self._reference_to_num:
                # Show some example mappings
                examples = list(self._reference_to_num.items())[:5]
                print(f"  Example reference mappings: {dict(examples)}")
        else:
            print("Reference format: Not detected")
        print()
        
        # NCBI RefSeq mappings
        if self._nc_to_num:
            print(f"NCBI RefSeq mappings: Available (build {self.build})")
            if self._num_to_nc:
                # Show some example mappings
                examples = list(self._num_to_nc.items())[:5]
                print(f"  Example NC mappings: {dict(examples)}")
        else:
            print("NCBI RefSeq mappings: Not available (build not specified or not supported)")
        print()
        
        # Middle layer (numeric) mappings
        print(f"Middle layer (numeric) mappings: Available")
        if self._num_to_str:
            # Show some example mappings
            examples = {k: v for k, v in list(self._num_to_str.items())[:5]}
            print(f"  Example numeric→string: {examples}")
        print()
        
        # Example conversions if formats are detected
        if self._sumstats_format and self._reference_format:
            print("Example conversions:")
            test_chrs = [1, 2, 23] if self._max_chr >= 23 else [1, 2]
            for num in test_chrs:
                if num in self._num_to_sumstats and num in self._num_to_reference:
                    sumstats_chr = self._num_to_sumstats[num]
                    ref_chr = self._num_to_reference[num]
                    print(f"  {sumstats_chr} (sumstats) ↔ {num} ↔ {ref_chr} (reference)")
        print("=" * 70)
