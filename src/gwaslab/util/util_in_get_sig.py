from typing import TYPE_CHECKING, Union, Optional, List, Dict, Any, Tuple
from gwaslab.qc.qc_check_datatype import check_datatype
import pandas as pd
import numpy as np
import scipy as sp
import gc
from os import path
from gwaslab.info.g_Log import Log

from gwaslab.bd.bd_common_data import get_chr_to_number
from gwaslab.bd.bd_common_data import get_number_to_chr
from gwaslab.bd.bd_common_data import get_chr_to_NC
from gwaslab.bd.bd_chromosome_mapper import ChromosomeMapper
from gwaslab.io.io_gtf import gtf_to_protein_coding
from gwaslab.io.io_gtf import gtf_to_all_gene
from gwaslab.io.io_gtf import read_gtf
from gwaslab.bd.bd_download import check_and_download
from gwaslab.qc.qc_build import _process_build
from gwaslab.qc.qc_fix_sumstats import check_dataframe_shape
from gwaslab.qc.qc_build import _check_build
from gwaslab.util.util_in_correct_winnerscurse import wc_correct
from gwaslab.util.util_ex_gwascatalog import gwascatalog_trait
from gwaslab.extension.gwascatalog import GWASCatalogClient
import gwaslab as gl
from gwaslab.util.util_in_fill_data import fill_p
from gwaslab.util.util_in_get_density import _get_signal_density2
from gwaslab.qc.qc_decorator import with_logging

if TYPE_CHECKING:
    from gwaslab.g_Sumstats import Sumstats

@with_logging(
        start_to_msg="extract lead variants",
        finished_msg="extracting lead variants",
        start_cols=["CHR","POS"],
        start_function=".get_lead()",
        fix=True
)
def _get_sig(insumstats_or_dataframe: Union['Sumstats', pd.DataFrame],
           variant_id: str = "SNPID",
           chrom: str = "CHR",
           pos: str = "POS",
           p: str = "P",
           mlog10p: str = "MLOG10P",
           scaled: bool = False,
           use_p: bool = False,
           windowsizekb: int = 500,
           sig_level: float = 5e-8,
           log: Log = Log(),
           xymt: List[str] = ["X","Y","MT"],
           anno: bool = False,
           wc_correction: bool = False,
           build: str = "19",
           source: str = "ensembl",
           gtf_path: Optional[str] = None,
           verbose: bool = True) -> Optional[pd.DataFrame]:
    """
    Extract lead variants by P values using a sliding window approach with significance thresholding.

    This function identifies lead variants from summary statistics using a sliding window 
    algorithm based on either -log10(p-values) or p-values. It prioritizes -log10(p-values) 
    if available, otherwise falls back to p-values. It handles data preprocessing, 
    significance filtering, and optional gene annotation and Winner's Curse correction.

    Parameters
    ----------
    insumstats_or_dataframe : Sumstats or pd.DataFrame
        Sumstats object or DataFrame to process.
    windowsizekb : int, default=500
        Window size in kilobases for lead variant identification, default=500
    sig_level : float, default=5e-8
        Significance threshold for variant selection, default=5e-8
    xymt : list, default=["X","Y","MT"]
        List of non-autosomal chromosome identifiers
    wc_correction : bool, default=False
        If True, apply Winner's Curse correction to effect sizes

    Returns
    -------
    pandas.DataFrame
        DataFrame containing significant lead variants with:
        - Original summary statistics columns
        - Annotated gene names (if anno=True)
        - Winner's Curse corrected BETA values (if wc_correction=True)
        - Additional metadata columns

    Notes
    -----
    The function performs multiple steps:
    1. Data validation and preprocessing
    2. Significance filtering using specified threshold
    3. Sliding window lead variant selection
    4. Optional gene annotation using Ensembl/RefSeq
    5. Optional Winner's Curse correction

    When no significant variants are found, returns None after logging a message.
    """
    import pandas as pd
    # Handle both DataFrame and Sumstats object
    if isinstance(insumstats_or_dataframe, pd.DataFrame):
        insumstats = insumstats_or_dataframe
    else:
        insumstats = insumstats_or_dataframe.data

    log.write(" -Processing "+str(len(insumstats))+" variants...", verbose=verbose)
    log.write(" -Significance threshold :", sig_level, verbose=verbose)
    log.write(" -Sliding window size:", str(windowsizekb) ," kb", verbose=verbose)
    
    #if "SNPID" in insumstats.columns:
    #    id = "SNPID"
    #else:
    #    id = "rsID"
#
    ##load data
    #sumstats=insumstats.loc[~insumstats[id].isna(),:].copy()
    
    sumstats = insumstats.copy()
    #convert chrom to int
    if sumstats[chrom].dtype in ["object",str,pd.StringDtype]:
        chr_to_num = get_chr_to_number(out_chr=True,xymt=["X","Y","MT"])
        sumstats[chrom]=sumstats[chrom].map(chr_to_num)
    
    # make sure the dtype is integer
    sumstats[chrom] = np.floor(pd.to_numeric(sumstats[chrom], errors='coerce')).astype('Int64')
    sumstats[pos] = np.floor(pd.to_numeric(sumstats[pos], errors='coerce')).astype('Int64')
    
    #create internal uniqid
    sumstats["__ID"] = range(len(sumstats))
    id_col = "__ID"

    #extract all significant variants
    ## use mlog10p first
    if use_p==False and (mlog10p in sumstats.columns):
        log.write(" -Using {} for extracting lead variants...".format(mlog10p),verbose=verbose)
        sumstats_sig = sumstats.loc[sumstats[mlog10p]>= -np.log10(sig_level),:].copy()
        sumstats_sig.loc[:,"__SCALEDP"] = -pd.to_numeric(sumstats_sig[mlog10p], errors='coerce')
    else:
        #use P
        log.write(" -Using {} for extracting lead variants...".format(p),verbose=verbose)         
        sumstats[p] = pd.to_numeric(sumstats[p], errors='coerce')
        sumstats_sig = sumstats.loc[sumstats[p]<=sig_level,:].copy()
        sumstats_sig.loc[:,"__SCALEDP"] = pd.to_numeric(sumstats_sig[p], errors='coerce')
    log.write(" -Found "+str(len(sumstats_sig))+" significant variants in total...", verbose=verbose)

    sumstats_sig = sumstats_sig.sort_values([chrom,pos])
    if len(sumstats_sig) == 0:
        log.write(" -No lead snps at given significance threshold!", verbose=verbose)
        return pd.DataFrame(columns=sumstats_sig.columns)

    p="__SCALEDP"
    sig_index_list = _collect_leads_generic(sumstats_sig, chrom, pos, windowsizekb, p, id_col, maximize=False)

    log.write(" -Identified "+str(len(sig_index_list))+" lead variants!", verbose=verbose)
    
    # drop internal __SCALEDP
    sumstats_sig = sumstats_sig.drop("__SCALEDP",axis=1)
    
    # extract the lead variants
    # Note: sig_index_list contains __ID values (internal IDs), not the original variant_id column values
    # So we need to filter using __ID, not the variant_id parameter
    output = sumstats_sig.loc[sumstats_sig[id_col].isin(sig_index_list),:].copy()

    # annotate GENENAME
    if anno is True and len(output)>0:
        log.write(" -Annotating variants using references:{}".format(source), verbose=verbose)
        log.write(" -Annotating variants using references based on genome build:{}".format(build), verbose=verbose)
        
        output = _anno_gene(
               output,
               id=variant_id,
               chrom=chrom,
               pos=pos,
               log=log,
               build=build,
               source=source,
               gtf_path=gtf_path,
               verbose=verbose)
        
    # drop internal id
    output = output.drop("__ID",axis=1)
    if wc_correction==True:
        log.write(" -Conducting Winner's Curse correction for BETA...", verbose=verbose)
        log.write(" -Referece: Zhong, H., & Prentice, R. L. (2008). Bias-reduced estimators and confidence intervals for odds ratios in genome-wide association studies. Biostatistics, 9(4), 621-634.")
        output["BETA_WC"] = output[["BETA","SE"]].apply(lambda x: wc_correct(x["BETA"],x["SE"],sig_level),axis=1)

    return output.copy()


@with_logging(
        start_to_msg="extract top variants by metric",
        finished_msg="extracting top variants",
        start_cols=["CHR","POS"],
        start_function=".get_top()",
        fix=True
)
def _get_top(
    insumstats_or_dataframe: Union['Sumstats', pd.DataFrame],
    variant_id: str = "SNPID",
    chrom: str = "CHR",
    pos: str = "POS",
    by: str = "DENSITY",
    threshold: Optional[float] = None,
    windowsizekb: int = 500,
    bwindowsizekb: int = 100,
    log: Log = Log(),
    xymt: List[str] = ["X","Y","MT"],
    anno: bool = False,
    build: str = "19",
    source: str = "ensembl",
    gtf_path: Optional[str] = None,
    verbose: bool = True
) -> Optional[pd.DataFrame]:
    """
    Extract top variants by maximizing a metric within sliding windows. (used for get top density variants)

    This function identifies top variants by selecting, within each
    contiguous window on a chromosome, the variant with the highest value
    of a specified column (e.g., `DENSITY`). It follows the same windowing
    logic as `getsig`, but does not rely on `P` or `MLOG10P`.

    Parameters
    ----------
    insumstats_or_dataframe : Sumstats or pd.DataFrame
        Sumstats object or DataFrame to process.
    by : str, default="DENSITY"
        Column name whose values are maximized to choose leads.
    threshold : float or None, default=None
        If provided, only variants with `by` >= `threshold` are considered. Deafult threshold is the median of maximum values of each chormosome.
    windowsizekb : int,
        Sliding window size in kilobases used to determine locus boundaries.  default=500
    bwindowsizekb : int, 
        Window size for calculating density. default=100
    anno : bool, default=False
        If True, annotate output with nearest gene names.

    Returns
    -------
    pandas.DataFrame or None
        DataFrame containing the selected lead variants. Returns None if
        no variants have valid values in the `by` column.
    """
    import pandas as pd
    # Handle both DataFrame and Sumstats object
    if isinstance(insumstats_or_dataframe, pd.DataFrame):
        insumstats = insumstats_or_dataframe
    else:
        insumstats = insumstats_or_dataframe.data
    log.write(" -Sliding window size for extracting top:", str(windowsizekb), " kb", verbose=verbose)
    sumstats = insumstats.copy()
    if sumstats[chrom].dtype in ["object",str,pd.StringDtype]:
        chr_to_num = get_chr_to_number(out_chr=True,xymt=["X","Y","MT"])
        sumstats[chrom]=sumstats[chrom].map(chr_to_num)
    sumstats[chrom] = np.floor(pd.to_numeric(sumstats[chrom], errors='coerce')).astype('Int64')
    sumstats[pos] = np.floor(pd.to_numeric(sumstats[pos], errors='coerce')).astype('Int64')

    if (by == "DENSITY") and ("DENSITY" not in sumstats.columns):
        sumstats = _get_signal_density2(insumstats=sumstats,
                                     snpid=variant_id,
                                     chrom=chrom,
                                     pos=pos,
                                     bwindowsizekb=bwindowsizekb,
                                     log=log,
                                     verbose=verbose)
    if by not in sumstats.columns:
        raise ValueError("Please make sure {} column is in input sumstats.".format(by))

    sumstats["__ID"] = range(len(sumstats))
    id_col = "__ID"
    sumstats_sig = sumstats.loc[~sumstats[by].isna(),:].copy()
    sumstats_sig["__BYNUM"] = pd.to_numeric(sumstats_sig[by], errors='coerce')
    if threshold is None:
        chr_max = sumstats_sig.groupby(chrom, sort=False)["__BYNUM"].max()
        threshold = float(chr_max.median()) if len(chr_max) > 0 else None
    if threshold is not None:
        sumstats_sig = sumstats_sig.loc[sumstats_sig["__BYNUM"] >= threshold, :]
        log.write(" -Using {} threshold: {}".format(by, threshold), verbose=verbose)
    sumstats_sig = sumstats_sig.sort_values([chrom,pos])
    if len(sumstats_sig) == 0:
        log.write(" -No variants passing {} threshold or valid metric".format(by), verbose=verbose)
        return None

    sig_index_list = _collect_leads_generic(sumstats_sig, chrom, pos, windowsizekb, "__BYNUM", id_col, maximize=True)
    log.write(" -Identified "+str(len(sig_index_list))+" top variants!", verbose=verbose)

    sumstats_sig = sumstats_sig.drop(["__BYNUM"],axis=1)
    output = sumstats_sig.loc[sumstats_sig[id_col].isin(sig_index_list),:].copy()

    if anno is True and len(output)>0:
        output = _anno_gene(output,
               id=variant_id,
               chrom=chrom,
               pos=pos,
               log=log,
               build=build,
               source=source,
               gtf_path=gtf_path,
               verbose=verbose)
    output = output.drop("__ID",axis=1)
    return output.copy()

def _collect_leads_generic(
    df: pd.DataFrame, 
    chrom: str, 
    pos: str, 
    windowsizekb: int, 
    score_col: str, 
    id_col: str, 
    maximize: bool = True
) -> List[Any]:
    # Return empty immediately if the input is empty
    if len(df) == 0:
        return []
    # Convert window size from kb to base pairs
    threshold = windowsizekb * 1000
    # Identify starts of new clusters when chromosome changes
    # Fill NA with True because first row (where shift() gives NaN) should start a new cluster
    new_chr = df[chrom].ne(df[chrom].shift()).fillna(True)
    # Compute position difference within chromosome; reset at new chromosome
    pos_diff = df[pos].diff().mask(new_chr, 0).fillna(0)
    # A new cluster starts on chromosome change or if gap exceeds threshold
    breaks = new_chr | (pos_diff > threshold)
    # Assign a cluster id by cumulative sum of break flags
    cluster_id = breaks.cumsum()
    # Choose the lead index per cluster by maximizing or minimizing the score
    # Check if score_col has any valid (non-NaN) values
    if df[score_col].notna().sum() == 0:
        return []
    # Add cluster_id as a temporary column to ensure proper alignment with DataFrame index
    df['__CLUSTER_ID'] = cluster_id
    try:
        if maximize:
            idx = df.groupby('__CLUSTER_ID', sort=False)[score_col].idxmax()
        else:
            idx = df.groupby('__CLUSTER_ID', sort=False)[score_col].idxmin()
        # Return the corresponding id values for selected indices
        # Use idx.values (row indices) not idx.index (cluster_ids) for lookup
        # Filter out any NaN values that might occur if score_col has all NaN in a cluster
        valid_idx = idx.dropna()
        if len(valid_idx) == 0:
            return []
        result = df.loc[valid_idx.values, id_col].tolist()
    finally:
        # Clean up temporary column
        df.drop('__CLUSTER_ID', axis=1, inplace=True)
    return result


class SimpleGenome:
    """
    Simple replacement for pyensembl.Genome that uses local GTF files.
    Provides gene_names_at_locus and gene_ids_at_locus methods.
    """
    def __init__(self, gtf_path_or_url, reference_name=None, annotation_name=None):
        """
        Initialize with a GTF file path.
        
        Parameters
        ----------
        gtf_path_or_url : str
            Path to GTF file (can be gzipped)
        reference_name : str, optional
            Reference name (for compatibility, not used)
        annotation_name : str, optional
            Annotation name (for compatibility, not used)
        """
        self.gtf_path = gtf_path_or_url
        self._genes_df = None
        self._indexed = False
    
    def index(self):
        """Load and index the GTF file. This is called automatically on first query."""
        if not self._indexed:
            # Load GTF file and filter to gene features
            self._genes_df = read_gtf(
                self.gtf_path,
                features={'gene'},
                usecols=['seqname', 'start', 'end', 'gene_id', 'gene_name'],
                expand_attribute_column=True
            )
            # Ensure seqname is string for consistent matching
            self._genes_df['seqname'] = self._genes_df['seqname'].astype(str)
            self._indexed = True
    
    def _ensure_indexed(self):
        """Ensure the GTF is loaded and indexed."""
        if not self._indexed:
            self.index()
    
    def _normalize_contig(self, contig):
        """Normalize chromosome/contig name for matching."""
        # Use ChromosomeMapper for normalization
        # Normalize to string format (preserves NC format, converts numeric to string)
        # Default to human species if not specified
        mapper = ChromosomeMapper(species="homo sapiens", verbose=False)
        # Detect format and convert to string
        mapper.detect_sumstats_format(pd.Series([contig]))
        normalized = mapper.to_string(contig)
        return str(normalized)
    
    def gene_names_at_locus(self, contig, position):
        """
        Get gene names at a specific genomic position.
        
        Parameters
        ----------
        contig : str
            Chromosome/contig name (e.g., "1", "X", "chr1")
        position : int
            Genomic position (1-based)
        
        Returns
        -------
        list
            List of gene names at the position (empty if none)
        """
        self._ensure_indexed()
        contig_norm = self._normalize_contig(contig)
        
        # Filter genes that overlap this position
        mask = (
            (self._genes_df['seqname'] == contig_norm) &
            (self._genes_df['start'] <= position) &
            (self._genes_df['end'] >= position)
        )
        matching_genes = self._genes_df[mask]
        
        if len(matching_genes) == 0:
            return []
        
        # Get unique gene names, filtering out None/NaN
        gene_names = matching_genes['gene_name'].dropna().unique().tolist()
        # Filter out empty strings
        gene_names = [g for g in gene_names if g and str(g).strip()]
        return gene_names
    
    def gene_ids_at_locus(self, contig, position):
        """
        Get gene IDs at a specific genomic position.
        
        Parameters
        ----------
        contig : str
            Chromosome/contig name (e.g., "1", "X", "chr1")
        position : int
            Genomic position (1-based)
        
        Returns
        -------
        list
            List of gene IDs at the position (empty if none)
        """
        self._ensure_indexed()
        contig_norm = self._normalize_contig(contig)
        
        # Filter genes that overlap this position
        mask = (
            (self._genes_df['seqname'] == contig_norm) &
            (self._genes_df['start'] <= position) &
            (self._genes_df['end'] >= position)
        )
        matching_genes = self._genes_df[mask]
        
        if len(matching_genes) == 0:
            return []
        
        # Get unique gene IDs, filtering out None/NaN
        gene_ids = matching_genes['gene_id'].dropna().unique().tolist()
        # Filter out empty strings
        gene_ids = [g for g in gene_ids if g and str(g).strip()]
        return gene_ids
    
    def closest_gene(self, chrom, pos, use_gene_id=False, max_distance=1000000):
        """
        Find the closest gene to a given chromosome and position.
        
        Parameters
        ----------
        chrom : str or int
            Chromosome number or name (e.g., "1", "X", 23, or NC_ format for RefSeq)
        pos : int
            Genomic position (1-based)
        use_gene_id : bool, default False
            If True, use gene_id instead of gene_name (for RefSeq)
        max_distance : int, default 1000000
            Maximum distance to search for genes (1Mb default)
        
        Returns
        -------
        tuple
            (distance, gene_names) where:
            - distance: int, distance to gene (0 if within gene, negative if upstream, positive if downstream)
            - gene_names: str, comma-separated gene names (or "intergenic" if no gene found)
        """
        self._ensure_indexed()
        
        # Normalize chromosome name
        contig_norm = self._normalize_contig(str(chrom))
        
        # Filter genes on this chromosome
        chr_genes = self._genes_df[self._genes_df['seqname'] == contig_norm].copy()
        
        if len(chr_genes) == 0:
            return max_distance, "intergenic"
        
        # Check if position is within any gene
        within_genes = chr_genes[
            (chr_genes['start'] <= pos) & (chr_genes['end'] >= pos)
        ]
        
        if len(within_genes) > 0:
            # Position is within gene(s)
            if use_gene_id:
                gene_col = 'gene_id'
            else:
                gene_col = 'gene_name'
            gene_names = within_genes[gene_col].dropna().unique().tolist()
            gene_names = [str(g).strip() for g in gene_names if g and str(g).strip()]
            if gene_names:
                return 0, ",".join(gene_names)
        
        # Calculate distances to all genes
        # Distance is negative if upstream, positive if downstream
        # If position is before gene start: distance = start - pos (positive = downstream)
        # If position is after gene end: distance = pos - end (positive = downstream)
        # We want the minimum absolute distance
        distances = []
        for _, gene in chr_genes.iterrows():
            if pos < gene['start']:
                # Position is upstream of gene
                dist = gene['start'] - pos
            elif pos > gene['end']:
                # Position is downstream of gene
                dist = pos - gene['end']
            else:
                # Shouldn't happen (we checked above), but handle it
                dist = 0
            
            if dist <= max_distance:
                distances.append((dist, gene))
        
        if not distances:
            return max_distance, "intergenic"
        
        # Find the closest gene(s)
        min_dist = min(d[0] for d in distances)
        closest_genes = [d[1] for d in distances if d[0] == min_dist]
        
        # Get gene names
        if use_gene_id:
            gene_col = 'gene_id'
        else:
            gene_col = 'gene_name'
        
        gene_names = [str(g[gene_col]).strip() for g in closest_genes 
                     if g[gene_col] and str(g[gene_col]).strip()]
        
        if not gene_names:
            return min_dist, "intergenic"
        
        # Determine if upstream or downstream
        # Check the first closest gene to determine direction
        first_gene = closest_genes[0]
        if pos < first_gene['start']:
            # Upstream
            distance = -min_dist
        else:
            # Downstream
            distance = min_dist
        
        return distance, ",".join(gene_names)
    
    def closest_genes_vectorized(self, chrom, positions, use_gene_id=False, max_distance=1000000):
        """
        Find closest genes for multiple positions on the same chromosome (vectorized).
        
        Parameters
        ----------
        chrom : str or int
            Chromosome number or name
        positions : array-like
            Array of genomic positions (1-based)
        use_gene_id : bool, default False
            If True, use gene_id instead of gene_name
        max_distance : int, default 1000000
            Maximum distance to search for genes
        
        Returns
        -------
        pandas.DataFrame
            DataFrame with columns ['distance', 'gene_names'] and same index as positions
        """
        self._ensure_indexed()
        
        # Normalize chromosome name
        contig_norm = self._normalize_contig(str(chrom))
        
        # Filter genes on this chromosome
        chr_genes = self._genes_df[self._genes_df['seqname'] == contig_norm].copy()
        
        if len(chr_genes) == 0:
            # No genes on this chromosome
            result = pd.DataFrame({
                'distance': [max_distance] * len(positions),
                'gene_names': ['intergenic'] * len(positions)
            }, index=range(len(positions)))
            return result
        
        # Convert positions to numpy array for vectorized operations
        positions = np.asarray(positions, dtype=np.int64)
        n_positions = len(positions)
        
        # Initialize results
        distances = np.full(n_positions, max_distance, dtype=np.int64)
        gene_names_list = ['intergenic'] * n_positions
        
        # Get gene column name
        if use_gene_id:
            gene_col = 'gene_id'
        else:
            gene_col = 'gene_name'
        
        # Convert gene coordinates to numpy arrays for vectorized operations
        gene_starts = chr_genes['start'].values
        gene_ends = chr_genes['end'].values
        gene_names_arr = chr_genes[gene_col].fillna('').astype(str).str.strip().values
        
        # Vectorized check for positions within genes
        # Create a 2D array: positions (rows) x genes (columns)
        # positions[:, None] creates a column vector, gene_starts/gene_ends are row vectors
        positions_2d = positions[:, None]
        within_mask = (positions_2d >= gene_starts) & (positions_2d <= gene_ends)
        
        # For each position, find genes it overlaps
        for pos_idx in range(n_positions):
            overlapping_genes = np.where(within_mask[pos_idx])[0]
            
            if len(overlapping_genes) > 0:
                # Position is within gene(s)
                distances[pos_idx] = 0
                names = [gene_names_arr[i] for i in overlapping_genes if gene_names_arr[i]]
                if names:
                    gene_names_list[pos_idx] = ",".join(names)
        
        # For positions not within genes, find closest gene using vectorized operations
        not_within_mask = distances != 0
        not_within_indices = np.where(not_within_mask)[0]
        
        if len(not_within_indices) > 0:
            not_within_positions = positions[not_within_indices]
            
            # For each position, calculate distances to all genes
            for pos_idx, pos in zip(not_within_indices, not_within_positions):
                # Calculate distances: upstream (pos < start) or downstream (pos > end)
                upstream_dist = np.where(pos < gene_starts, gene_starts - pos, np.inf)
                downstream_dist = np.where(pos > gene_ends, pos - gene_ends, np.inf)
                
                # Minimum distance for each gene
                gene_distances = np.minimum(upstream_dist, downstream_dist)
                
                # Find minimum distance
                min_dist = np.min(gene_distances)
                
                if min_dist < max_distance:
                    # Find all genes at minimum distance
                    closest_gene_indices = np.where(gene_distances == min_dist)[0]
                    
                    # Get gene names
                    names = [gene_names_arr[i] for i in closest_gene_indices if gene_names_arr[i]]
                    
                    if names:
                        # Determine direction (check first closest gene)
                        first_idx = closest_gene_indices[0]
                        if pos < gene_starts[first_idx]:
                            distances[pos_idx] = -int(min_dist)
                        else:
                            distances[pos_idx] = int(min_dist)
                        gene_names_list[pos_idx] = ",".join(names)
        
        # Create result DataFrame
        result = pd.DataFrame({
            'distance': distances,
            'gene_names': gene_names_list
        }, index=range(n_positions))
        
        return result


Genome = SimpleGenome


def closest_gene(
    x: pd.Series, 
    data: Any, 
    chrom: str = "CHR", 
    pos: str = "POS", 
    source: str = "ensembl", 
    build: str = "19"
) -> Tuple[int, str]:
    """
    Find the closest gene to a variant position.
    
    Parameters
    ----------
    x : pandas.Series
        Row from DataFrame containing variant information
    data : SimpleGenome
        SimpleGenome instance with loaded GTF data
    chrom : str, default "CHR"
        Column name for chromosome
    pos : str, default "POS"
        Column name for position
    source : str, default "ensembl"
        Data source ("ensembl" or "refseq")
    build : str, default "19"
        Genome build version
    
    Returns
    -------
    tuple
        (distance, gene_names) where distance is int and gene_names is str
    """
    use_gene_id = (source == "refseq")
    
    # Handle chromosome conversion for RefSeq
    if source == "refseq":
        chrom_str = get_number_to_chr()[x[chrom]]
        chrom_nc = get_chr_to_NC(build=build)[chrom_str]
        return data.closest_gene(chrom_nc, int(x[pos]), use_gene_id=use_gene_id)
    else:
        return data.closest_gene(x[chrom], int(x[pos]), use_gene_id=use_gene_id)


@with_logging(
        start_to_msg="annotate variants with nearest gene name(s)",
        finished_msg="annotating variants with nearest gene name(s) successfully!"
)
def _anno_gene(
    insumstats_or_dataframe: Union['Sumstats', pd.DataFrame],
    id: str = "SNPID",
    chrom: str = "CHR",
    pos: str = "POS",
    log: Log = Log(),
    build: str = "19",
    source: str = "ensembl",
    gtf_path: Optional[str] = None,
    verbose: bool = True
) -> pd.DataFrame:
    import pandas as pd
    # Handle both DataFrame and Sumstats object
    if isinstance(insumstats_or_dataframe, pd.DataFrame):
        insumstats = insumstats_or_dataframe
    else:
        insumstats = insumstats_or_dataframe.data
    
    # Auto-detect ID column if using default and column doesn't exist
    if id == "SNPID" and id not in insumstats.columns:
        if "rsID" in insumstats.columns:
            id = "rsID"
    
    build = _process_build(build, log=log,verbose=verbose)
    output = insumstats.copy()
    
    # Determine GTF path and prepare genome object
    if source == "ensembl":
        if build=="19":
            log.write(" -Assigning Gene name using ensembl_hg19_gtf for protein coding genes", verbose=verbose)
            if gtf_path is None:
                gtf_path = check_and_download("ensembl_hg19_gtf")
                gtf_path = gtf_to_protein_coding(gtf_path,log=log,verbose=verbose)
            else:
                log.write(" -Using user-provided gtf:{}".format(gtf_path))
                gtf_path = gtf_to_all_gene(gtf_path,log=log,verbose=verbose)
            reference_name = 'GRCh37'
        elif build=="38":
            log.write(" -Assigning Gene name using ensembl_hg38_gtf for protein coding genes", verbose=verbose)
            if gtf_path is None:
                gtf_path = check_and_download("ensembl_hg38_gtf")
                gtf_path = gtf_to_protein_coding(gtf_path,log=log,verbose=verbose)
            else:
                log.write(" -Using user-provided gtf:{}".format(gtf_path))
                gtf_path = gtf_to_all_gene(gtf_path,log=log,verbose=verbose)
            reference_name = 'GRCh38'
        else:
            raise ValueError(f"Unsupported build: {build}")
        use_gene_id = False
    elif source == "refseq":
        if build=="19":
            log.write(" -Assigning Gene name using NCBI refseq latest GRCh37 for protein coding genes", verbose=verbose)
            if gtf_path is None:
                gtf_path = check_and_download("refseq_hg19_gtf")
                gtf_path = gtf_to_protein_coding(gtf_path,log=log,verbose=verbose)
            else:
                log.write(" -Using user-provided gtf:{}".format(gtf_path))
                gtf_path = gtf_to_all_gene(gtf_path,log=log,verbose=verbose)
            reference_name = 'GRCh37'
        elif build=="38":
            log.write(" -Assigning Gene name using NCBI refseq latest GRCh38 for protein coding genes", verbose=verbose)
            if gtf_path is None:
                gtf_path = check_and_download("refseq_hg38_gtf")
                gtf_path = gtf_to_protein_coding(gtf_path,log=log,verbose=verbose)
            else:
                log.write(" -Using user-provided gtf:{}".format(gtf_path))
                gtf_path = gtf_to_all_gene(gtf_path,log=log,verbose=verbose)
            reference_name = 'GRCh38'
        else:
            raise ValueError(f"Unsupported build: {build}")
        use_gene_id = True
    else:
        raise ValueError(f"Unsupported source: {source}")
    
    # Check if input is empty
    if len(output) == 0:
        output["LOCATION"] = ""
        output["GENE"] = ""
        return output
    
    # Initialize genome object
    data = Genome(
        reference_name=reference_name,
        annotation_name=source.capitalize(),
        gtf_path_or_url=gtf_path)
    
    # Pre-index the genome (load GTF data once)
    log.write(" -Loading and indexing GTF file...", verbose=verbose)
    data.index()
    
    # Initialize result columns
    output["LOCATION"] = 0
    output["GENE"] = "Unknown"
    
    # Process by chromosome for efficiency
    log.write(" -Processing variants by chromosome...", verbose=verbose)
    
    # Get unique chromosomes
    unique_chroms = output[chrom].unique()
    
    # Process each chromosome
    for chrom_val in unique_chroms:
        # Get variants on this chromosome
        chrom_mask = output[chrom] == chrom_val
        chrom_indices = output[chrom_mask].index
        chrom_variants = output.loc[chrom_mask]
        
        if len(chrom_variants) == 0:
            continue
        
        # Convert chromosome name for query
        if source == "refseq":
            chrom_str = get_number_to_chr()[chrom_val]
            query_chrom = get_chr_to_NC(build=build)[chrom_str]
        else:
            query_chrom = chrom_val
        
        # Get positions
        positions = chrom_variants[pos].values
        
        # Use vectorized method to find closest genes
        results = data.closest_genes_vectorized(
            query_chrom, 
            positions, 
            use_gene_id=use_gene_id
        )
        
        # Update output using the original indices
        output.loc[chrom_indices, "LOCATION"] = results["distance"].values
        output.loc[chrom_indices, "GENE"] = results["gene_names"].values
        
        if verbose and len(unique_chroms) > 1:
            log.write(f"   -Processed {len(chrom_variants)} variants on chromosome {chrom_val}", verbose=verbose)
    
    # Replace empty strings and "intergenic" with "Unknown"
    output["GENE"] = output["GENE"].replace("","Unknown")
    output["GENE"] = output["GENE"].replace("intergenic","Unknown")
    
    return output

def _get_known_variants_from_gwascatalog(
    client: Any, 
    efo: str, 
    sig_level: float = 5e-8, 
    verbose: bool = True, 
    log: Log = Log()
) -> pd.DataFrame:
    """
    Retrieve known variants from GWAS Catalog API v2 for a given EFO trait.
    
    This function delegates to the GWASCatalogClient.get_known_variants_for_trait method,
    which handles all the GWAS Catalog API logic, then processes the results.
    
    Parameters
    ----------
    client : GWASCatalogClient
        GWAS Catalog API client instance
    efo : str
        EFO trait ID (e.g., "EFO_0001360"), MONDO ID (e.g., "MONDO_0005148"),
        or trait name (e.g., "duodenal ulcer")
    sig_level : float
        P-value threshold for filtering associations
    verbose : bool
        Whether to print log messages
    log : Log
        Logging object
        
    Returns
    -------
    pd.DataFrame
        DataFrame with columns: SNPID, CHR, POS, and other association metadata
    """
    
    # Delegate to the client method which handles all GWAS Catalog API logic
    knownsig = client.get_known_variants_for_trait(efo=efo, sig_level=sig_level, verbose=verbose)
    
    if len(knownsig) == 0:
        return pd.DataFrame()
    
    # Remove duplicates based on CHR:POS
    before_dedup = len(knownsig)
    knownsig = knownsig.drop_duplicates(subset=['CHR', 'POS'], keep='first')
    if before_dedup > len(knownsig):
        log.write(" -Removed {} duplicate variants (kept {} unique)".format(before_dedup - len(knownsig), len(knownsig)), verbose=verbose)
    
    # Create Sumstats object for consistency with old format and process coordinates
    known_Sumstats = gl.Sumstats(knownsig.copy(), fmt="gwaslab", 
                                 other=['REPORT_GENENAME', 'TRAIT', 'STUDY', 'PUBMEDID', 'AUTHOR'],
                                 verbose=False)
    known_Sumstats.fix_pos(verbose=False)
    known_Sumstats.fix_chr(verbose=False)
    known_Sumstats.sort_coordinate(verbose=False)
    
    log.write(" -Processed {} unique variants from GWAS Catalog".format(len(known_Sumstats.data)), verbose=verbose)
    
    return known_Sumstats.data

@with_logging(
        start_to_msg="check if lead variants are known",
        finished_msg="checking if lead variants are known",
        start_cols=["CHR","POS"],
        start_function=".get_novel()"
)
def _get_novel(
    insumstats_or_dataframe: Union['Sumstats', pd.DataFrame],
    variant_id: str = "SNPID",
    chrom: str = "CHR",
    pos: str = "POS",
    p: str = "P",
    use_p: bool = False,
    known: Union[bool, pd.DataFrame, str] = False,
    efo: Union[bool, str, List[str]] = False,
    only_novel: bool = False,
    group_key: Optional[str] = None,
    if_get_lead: bool = True,
    windowsizekb_for_novel: int = 1000,
    windowsizekb: int = 500,
    sig_level: float = 5e-8,
    log: Log = Log(),
    xymt: List[str] = ["X","Y","MT"],
    anno: bool = False,
    wc_correction: bool = False,
    use_cache: bool = True,
    cache_dir: str = "./",
    build: str = "19",
    source: str = "ensembl",
    gwascatalog_source: str = "NCBI",
    output_known: bool = False,
    verbose: bool = True
) -> pd.DataFrame:
    """
    Identify novel variants by comparing against known variant databases.

    This function determines whether variants in summary statistics are novel by comparing them 
    against known variants from GWAS catalog or user-provided reference data. It handles 
    coordinate conversion, distance calculations, and group-based comparisons.

    Parameters
    ----------
    insumstats_or_dataframe : Sumstats or pd.DataFrame
        Sumstats object or DataFrame to process.
    known : pd.DataFrame or str, optional
        DataFrame or path to file containing known variants with CHR/POS columns
    efo : str or list, optional
        EFO trait(s) for querying GWAS catalog
    only_novel : bool, default=False
        If True, return only novel variants
    group_key : str, optional
        Column name for grouping variants (e.g., trait/phenotype ID)
    if_get_lead : bool, default=True
        If True, first extract lead variants using getsig
    windowsizekb : int, default=500
        Window size (kb) for lead variant identification
    windowsizekb_for_novel : int, default=1000
        Distance threshold (kb) to define novelty

    Returns
    -------
    pandas.DataFrame or tuple
        If only_novel=False and output_known=False: DataFrame with all variants and NOVEL column
        If only_novel=True and output_known=False: DataFrame with only novel variants
        If output_known=True: tuple of (variants DataFrame, known variants DataFrame)
        The returned DataFrame includes a "NOVEL" column indicating novelty status.

    Notes
    -----
    The function performs multiple steps:
    1. Data validation and preprocessing
    2. Retrieval of known variants from GWAS catalog or user input
    3. Coordinate conversion and helper column creation (TCHR+POS)
    4. Distance calculations between variants
    5. Novelty determination based on distance threshold
    6. Grouped comparisons when group_key is provided

    When no significant variants are found, returns None after logging a message.
    """
    import pandas as pd
    # Handle both DataFrame and Sumstats object
    if isinstance(insumstats_or_dataframe, pd.DataFrame):
        insumstats = insumstats_or_dataframe
    else:
        insumstats = insumstats_or_dataframe.data

    if "SNPID" in insumstats.columns:
        variant_id = "SNPID"
    else:
        variant_id = "rsID"

    if if_get_lead == True:
        allsig = _get_sig(insumstats=insumstats,
            variant_id=variant_id,chrom=chrom,pos=pos,p=p,use_p=use_p,windowsizekb=windowsizekb,sig_level=sig_level,log=log,
            xymt=xymt,anno=anno,build=build,wc_correction=wc_correction, source=source,verbose=verbose)
    else:
        allsig = insumstats.copy()

    ############################################################################################
    knownsig = pd.DataFrame()
    if efo != False:
        # For GWAS catalog, checking if sumstats build is hg38
        _check_build(target_build="38" ,build=build ,log=log,verbose=verbose)
        
        # Use new GWAS Catalog API v2 client
        client = GWASCatalogClient(verbose=verbose, log=log)
        
        if type(efo) is not list:
            log.write("Start to retrieve data using EFO: {}...".format(efo), verbose=verbose)
            knownsig = _get_known_variants_from_gwascatalog(
                client=client,
                efo=efo,
                sig_level=sig_level,
                verbose=verbose,
                log=log
            )
        else:            
            knownsig=pd.DataFrame()
            log.write("Start to retrieve data using {} EFOs: {}...".format(len(efo),efo), verbose=verbose)

            for single_efo in efo:
                knownsig_single = _get_known_variants_from_gwascatalog(
                    client=client,
                    efo=single_efo,
                    sig_level=sig_level,
                    verbose=verbose,
                    log=log
                )
                if len(knownsig_single) > 0:
                    knownsig_single["EFOID"] = single_efo
                    knownsig = pd.concat([knownsig, knownsig_single],ignore_index=True)
        
        if len(knownsig) > 0:
            knownsig["CHR"] = knownsig["CHR"].astype("Int64")
            knownsig["POS"] = knownsig["POS"].astype("Int64")
            log.write(" -Retrieved {} associations from GWAS catalog.".format(len(knownsig)), verbose=verbose)
        else:
            log.write(" -No associations found in GWAS catalog for the specified EFO trait(s).", verbose=verbose)
    if type(known) is pd.DataFrame:
        knownsig_2 = known.copy()
        knownsig = pd.concat([knownsig, knownsig_2],ignore_index=True)
        knownsig["CHR"] = knownsig["CHR"].astype("Int64")
        knownsig["POS"] = knownsig["POS"].astype("Int64")
        if "SNPID" not in knownsig.columns:
            knownsig["SNPID"] =knownsig["CHR"].astype("string") + ":" + knownsig["POS"].astype("string")
    elif type(known) is str:
        knownsig_2 = pd.read_csv(known,sep="\s+",dtype={"CHR":"Int64","POS":"Int64"})
        knownsig = pd.concat([knownsig, knownsig_2],ignore_index=True)
        knownsig["CHR"] = knownsig["CHR"].astype("Int64")
        knownsig["POS"] = knownsig["POS"].astype("Int64")
        if "SNPID" not in knownsig.columns:
            knownsig["SNPID"] =knownsig["CHR"].astype("string") + ":" + knownsig["POS"].astype("string")
    
    if len(knownsig)<1:
        raise ValueError("Please input a dataframe of known loci or valid efo code")
    
    if group_key is not None:
        if (group_key not in allsig.columns) or (group_key not in knownsig.columns):
            raise ValueError("Please check if group_key is in both sumstats and list of known associations.")

    # create helper column TCHR+POS for knownsig and all sig
     ############################################################################################
    maxpos = insumstats["POS"].max()
    big_number = determine_big_number(maxpos)
    knownsig = add_tchr_pos(knownsig, chrom, pos, big_number)
    allsig = add_tchr_pos(allsig, chrom, pos, big_number)
    ############################################################################################   
    #sorting
    allsig = allsig.sort_values(by="TCHR+POS",ignore_index=True)
    knownsig = knownsig.sort_values(by="TCHR+POS",ignore_index=True)
    ############################################################################################
    if group_key is not None:
        number_of_groups_allsig = allsig[group_key].nunique()
        number_of_groups_known = knownsig[group_key].nunique()
        log.write(" -Number of groups in sumstats:{}".format(number_of_groups_allsig), verbose=verbose)
        log.write(" -Number of groups in reference:{}".format(number_of_groups_known), verbose=verbose)

    log.write(" -Lead variants in known loci:",len(knownsig), verbose=verbose)
    log.write(" -Checking the minimum distance between identified lead variants and provided known variants...", verbose=verbose)
    
    ############################################################################################
    if group_key is None:
        # get distance
        allsig = determine_distance(allsig, knownsig)
        # get other info 
        allsig = fill_meta_info_for_known(allsig, knownsig)
        ############################################################################################
        # determine if novel
        allsig = determine_novel(allsig, windowsizekb_for_novel)
        # determine location
        allsig = determine_location(allsig)
        # if not on same chromosome, distance set to pd.NA
        allsig = determine_if_same_chromosome(allsig, knownsig, maxpos)
        ############################################################################################
    else:
        #groups1 = set(allsig[group_key].unique())
        #groups2 = set(knownsig[group_key].unique())
        #common_group = groups1.intersection(groups2)
        
        #allsig_no_group = allsig.loc[~allsig[group_key].isin(common_group),:].copy()
        allsig_group = pd.DataFrame()

        for key in allsig[group_key].unique():
            allsig_single_group = allsig.loc[allsig[group_key]==key,:].copy()
            knownsig_single_group = knownsig.loc[knownsig[group_key]==key,:].copy()

            #if len(allsig_single_group) >0 and len(knownsig_single_group) >0:
            allsig_single_group = determine_distance(allsig_single_group, knownsig_single_group)
            # get other info 
            allsig_single_group = fill_meta_info_for_known(allsig_single_group, knownsig_single_group)
            
            # determine if novel
            allsig_single_group = determine_novel(allsig_single_group, windowsizekb_for_novel)
            
            # determine location
            allsig_single_group = determine_location(allsig_single_group)
            
            # if not on same chromosome, distance set to pd.NA
            allsig_single_group = determine_if_same_chromosome(allsig_single_group, knownsig_single_group, maxpos) 
            
            allsig_group = pd.concat([allsig_group, allsig_single_group], ignore_index=True)
        
        allsig = allsig_group
        #pd.concat([allsig_no_group, allsig_group], ignore_index=True)
        
    # drop helper column TCHR+POS
    allsig = allsig.drop(["TCHR+POS"], axis=1)
    
    try:
        allsig = allsig.where(~pd.isna(allsig), pd.NA)
    except:
        pass

    log.write(" -Identified ",len(allsig)-sum(allsig["NOVEL"])," known variants in current sumstats...", verbose=verbose)
    log.write(" -Identified ",sum(allsig["NOVEL"])," novel variants in current sumstats...", verbose=verbose)
    
    # how to return
    if only_novel is True:
        if output_known is True:
            return allsig.loc[allsig["NOVEL"],:], knownsig
        else:
            return allsig.loc[allsig["NOVEL"],:]
    else:
        if output_known is True:
            return allsig, knownsig
        else:
            return allsig
##################################################################################################################################################################################################

@with_logging(
        start_to_msg="check if variants are in cis or trans regions",
        finished_msg="checking if variants are in cis or trans regions",
        start_cols=["CHR","POS"],
        start_function=".check_cis()",
        must_kwargs=["group_key"]
)
def _check_cis(
    insumstats_or_dataframe: Union['Sumstats', pd.DataFrame],
    variant_id: Optional[str] = None,
    chrom: str = "CHR",
    pos: str = "POS",
    p: str = "P",
    use_p: bool = False,
    known: Union[bool, pd.DataFrame, str] = False,
    group_key: Optional[str] = None,
    if_get_lead: bool = False,
    windowsizekb: int = 500,
    sig_level: float = 5e-8,
    log: Log = Log(),
    xymt: List[str] = ["X","Y","MT"],
    anno: bool = False,
    build: str = "19",
    source: str = "ensembl",
    verbose: bool = True
) -> pd.DataFrame:
    import pandas as pd
    # Handle both DataFrame and Sumstats object
    if isinstance(insumstats_or_dataframe, pd.DataFrame):
        insumstats = insumstats_or_dataframe
    else:
        insumstats = insumstats_or_dataframe.data
    
    ##start function with col checking##########################################################
    
    # Auto-detect ID column if not provided
    if variant_id is None:
        if "SNPID" in insumstats.columns:
            variant_id = "SNPID"
        elif "rsID" in insumstats.columns:
            variant_id = "rsID"
        else:
            raise ValueError("Cannot find SNPID or rsID column in sumstats")
    
    if if_get_lead == True:
        allsig = _get_sig(insumstats=insumstats,
            variant_id=variant_id,chrom=chrom,pos=pos,p=p,use_p=use_p,windowsizekb=windowsizekb,sig_level=sig_level,log=log,
            xymt=xymt,anno=anno,build=build, source=source,verbose=verbose)
    else:
        allsig = insumstats.copy()

    ############################################################################################
    knownsig = pd.DataFrame()
    if type(known) is pd.DataFrame:
        knownsig_2 = known.copy()
        knownsig = pd.concat([knownsig, knownsig_2],ignore_index=True)
        knownsig["CHR"] = knownsig["CHR"].astype("Int64")
        knownsig["START"] = knownsig["START"].astype("Int64")
        knownsig["END"] = knownsig["END"].astype("Int64")
    elif type(known) is str:
        knownsig_2 = pd.read_csv(known,sep="\s+",dtype={"CHR":"Int64","POS":"Int64"})
        knownsig = pd.concat([knownsig, knownsig_2],ignore_index=True)
        knownsig["CHR"] = knownsig["CHR"].astype("Int64")
        knownsig["START"] = knownsig["START"].astype("Int64")
        knownsig["END"] = knownsig["END"].astype("Int64")
    
    if len(knownsig)<1:
        raise ValueError("Please input a dataframe of gene list with GENE, CHR, START, END.")
    
    if group_key is not None:
        if group_key not in knownsig.columns:
            raise ValueError("Please check if group_key is in both sumstats and list of known associations.")

    ############################################################################################
    if group_key is not None:
        number_of_groups_allsig = allsig[group_key].nunique()
        number_of_groups_known = knownsig[group_key].nunique()
        log.write(" -Number of groups in sumstats:{}".format(number_of_groups_allsig), verbose=verbose)
        log.write(" -Number of groups in reference:{}".format(number_of_groups_known), verbose=verbose)

    log.write(" -Checking if variants in cis/trans regions grouped by {}...".format(group_key), verbose=verbose)
    log.write(" -Window size in kb adding to start and end: {}...".format(windowsizekb), verbose=verbose)
    ############################################################################################
    #convert to  a dict
    reference_dict = {}
    for index,row in knownsig.iterrows():
        reference_dict[row[group_key]] = (row["CHR"], row["START"], row["END"] )
    ############################################################################################
    try:
        no_reference_avaialble = allsig.loc[~allsig[group_key].isin(reference_dict.keys()),group_key]
        if len(no_reference_avaialble)>0:
            log.write(" -Groups not in reference: {}".format( ",".join(no_reference_avaialble.unique())), verbose=verbose)
    except:
        pass

    #allsig["CIS/TRANS"] = allsig.apply(lambda x: determine_if_cis(x, group_key,windowsizekb, reference_dict), axis=1)
    cis_tuples = allsig.apply(lambda x: determine_if_cis2(x, group_key,windowsizekb, reference_dict), axis=1)
    allsig[["CIS/TRANS","REF_CHR","REF_START","REF_END"]] = pd.DataFrame(cis_tuples.tolist(), index=allsig.index)

    try:
        allsig = allsig.where(~pd.isna(allsig), pd.NA)
    except:
        pass
    
    try:
        number_of_cis = sum(allsig["CIS/TRANS"] == "Cis")
        number_of_trans = sum(allsig["CIS/TRANS"] == "Trans")
        number_of_noreference = sum(allsig["CIS/TRANS"] == "NoReference")
        log.write (" -Number of Cis variants: {}".format(number_of_cis),verbose=verbose)
        log.write (" -Number of Trans variants: {}".format(number_of_trans),verbose=verbose)
        log.write (" -Number of NoReference variants: {}".format(number_of_noreference),verbose=verbose)
    except:
        pass
    
    return allsig

###################################################################################################################################################################################################


def determine_big_number(maxpos, big_number = 1000000000):
    for i in range(7):
        if maxpos*10 >  big_number:
            big_number = int(big_number * 10)
        else:
            break
    return big_number

def add_tchr_pos(df, chrom, pos, big_number):
    df["TCHR+POS"]=df[chrom]*big_number + df[pos]
    return df

def fill_meta_info_for_known(allsig, knownsig):
    if len(allsig)==0 or len(knownsig)==0: return allsig
    if "SNPID" in knownsig.columns:
        knownids=knownsig["SNPID"].values
    if "PUBMEDID" in knownsig.columns:
        knownpubmedids=knownsig["PUBMEDID"].values
    if "AUTHOR" in knownsig.columns:
        knownauthor=knownsig["AUTHOR"].values
    if "EFOID" in knownsig.columns:
        knownefo=knownsig["EFOID"].values

    if "SNPID" in knownsig.columns:
        allsig["KNOWN_ID"] = allsig["TCHR+POS"].apply(lambda x:knownids[np.argmin(np.abs(knownsig["TCHR+POS"]-x))])    
    if "PUBMEDID" in knownsig.columns:
        allsig["KNOWN_PUBMED_ID"] = allsig["TCHR+POS"].apply(lambda x:knownpubmedids[np.argmin(np.abs(knownsig["TCHR+POS"]-x))])
    if "AUTHOR" in knownsig.columns:
        allsig["KNOWN_AUTHOR"] = allsig["TCHR+POS"].apply(lambda x:knownauthor[np.argmin(np.abs(knownsig["TCHR+POS"]-x))])
    if "EFOID" in knownsig.columns:
        allsig["KNOWN_EFOID"] = allsig["TCHR+POS"].apply(lambda x:knownefo[np.argmin(np.abs(knownsig["TCHR+POS"]-x))])    
    return allsig

def determine_if_cis(x, group_key,windowsizekb, reference_dict):
    if x[group_key] in reference_dict.keys():
        is_same_chr = str(reference_dict[x[group_key]][0]) == str(x["CHR"])
        is_large_than_start = int(reference_dict[x[group_key]][1]) - windowsizekb*1000 <= x["POS"]
        is_smaller_than_end = int(reference_dict[x[group_key]][2]) + windowsizekb*1000 >= x["POS"]               
        
        if  is_same_chr and is_large_than_start  and is_smaller_than_end:
            return "Cis"
        else: 
            return "Trans"
    else:
        return "NoReference"

def determine_if_cis2(x, group_key,windowsizekb, reference_dict):
    if x[group_key] in reference_dict.keys():
        is_same_chr = str(reference_dict[x[group_key]][0]) == str(x["CHR"])
        is_large_than_start = int(reference_dict[x[group_key]][1]) - windowsizekb*1000 <= x["POS"]
        is_smaller_than_end = int(reference_dict[x[group_key]][2]) + windowsizekb*1000 >= x["POS"]               
        
        if  is_same_chr and is_large_than_start  and is_smaller_than_end:
            return "Cis", int(reference_dict[x[group_key]][0]), int(reference_dict[x[group_key]][1]), int(reference_dict[x[group_key]][2])
        else: 
            return "Trans", int(reference_dict[x[group_key]][0]), int(reference_dict[x[group_key]][1]), int(reference_dict[x[group_key]][2])
    else:
        return "NoReference", pd.NA, pd.NA, pd.NA


def determine_distance(allsig, knownsig):
    if len(allsig)==0: 
        return allsig
    if len(knownsig)==0:
        allsig["DISTANCE_TO_KNOWN"] = pd.NA
        return allsig
    allsig["DISTANCE_TO_KNOWN"] = allsig["TCHR+POS"].apply(lambda x:min(knownsig["TCHR+POS"]-x, key=abs))
    return allsig

def determine_novel(allsig, windowsizekb_for_novel):
    if len(allsig)==0 or "DISTANCE_TO_KNOWN" not in allsig.columns:
        return allsig
    allsig["NOVEL"] = allsig["DISTANCE_TO_KNOWN"].abs() > windowsizekb_for_novel*1000
    allsig.loc[allsig["DISTANCE_TO_KNOWN"].isna(), "NOVEL"] = True
    return allsig

def determine_location(allsig):
    allsig["LOCATION_OF_KNOWN"]="NoReference"
    allsig.loc[ allsig["DISTANCE_TO_KNOWN"]== 0,"LOCATION_OF_KNOWN"] = "Same"
    allsig.loc[ allsig["DISTANCE_TO_KNOWN"] > 0 ,"LOCATION_OF_KNOWN"] = "Upstream"
    allsig.loc[ allsig["DISTANCE_TO_KNOWN"] < 0 ,"LOCATION_OF_KNOWN"] = "Downstream"
    return allsig

def determine_if_same_chromosome(allsig, knownsig, maxpos):
    if sum(allsig["DISTANCE_TO_KNOWN"].abs() > maxpos)>0:
        not_on_same_chromosome = allsig["DISTANCE_TO_KNOWN"].abs() > maxpos
        allsig.loc[ not_on_same_chromosome ,"DISTANCE_TO_KNOWN"] = pd.NA
        allsig.loc[ not_on_same_chromosome ,"LOCATION_OF_KNOWN"] = "NoneOnThisChr"
        if "SNPID" in knownsig.columns:
            allsig.loc[ not_on_same_chromosome ,"KNOWN_ID"] = pd.NA
        if "PUBMEDID" in knownsig.columns:
            allsig.loc[ not_on_same_chromosome ,"KNOWN_PUBMED_ID"] = pd.NA
        if "AUTHOR" in knownsig.columns:
            allsig.loc[ not_on_same_chromosome ,"KNOWN_AUTHOR"] = pd.NA
        if "EFOID" in knownsig.columns:
            allsig.loc[ not_on_same_chromosome ,"KNOWN_EFOID"] = pd.NA
    return allsig

@with_logging(
        start_to_msg="check if variant sets are overlapping with those in reference file",
        finished_msg="checking if variant sets are overlapping with those in reference file",
        start_cols=["CHR","POS"],
        start_function=".check_novel_set()",
        must_kwargs=["group_key"]
)
def _check_novel_set(insumstats_or_dataframe,
           variant_id=None,
           chrom="CHR",
           pos="POS",
           p="P",
           use_p=False,
           known=False,
           group_key=None,
           snpset="SNPSET",
           snpid="SNPID",
           if_get_lead = False,
           windowsizekb=500,
           sig_level=5e-8,
           log=Log(),
           xymt=["X","Y","MT"],
           anno=False,
           build="19",
           source="ensembl",
           verbose=True):
    import pandas as pd
    # Handle both DataFrame and Sumstats object
    if isinstance(insumstats_or_dataframe, pd.DataFrame):
        insumstats = insumstats_or_dataframe
    else:
        insumstats = insumstats_or_dataframe.data
    
    ##start function with col checking##########################################################
    
    # Auto-detect ID column if not provided
    if variant_id is None:
        if "SNPID" in insumstats.columns:
            variant_id = "SNPID"
        elif "rsID" in insumstats.columns:
            variant_id = "rsID"
        else:
            raise ValueError("Cannot find SNPID or rsID column in sumstats")
    
    # Auto-detect snpid for known set if using default and column doesn't exist
    if snpid == "SNPID" and snpid not in insumstats.columns:
        if "rsID" in insumstats.columns:
            snpid = "rsID"
    
    if if_get_lead == True:
        allsig = _get_sig(insumstats=insumstats,
            variant_id=variant_id,chrom=chrom,pos=pos,p=p,use_p=use_p,windowsizekb=windowsizekb,sig_level=sig_level,log=log,
            xymt=xymt,anno=anno,build=build, source=source,verbose=verbose)
    else:
        allsig = insumstats.copy()

    ############################################################################################
    knownsig = pd.DataFrame()
    if type(known) is pd.DataFrame:
        knownsig_2 = known.copy()
        knownsig = pd.concat([knownsig, knownsig_2],ignore_index=True)
        knownsig[snpid] = knownsig[snpid].astype("string")
        knownsig[snpset] = knownsig[snpset].astype("string")
        knownsig[group_key] = knownsig[group_key].astype("string")
    elif type(known) is str:
        knownsig_2 = pd.read_csv(known,sep="\s+",dtype={"CHR":"Int64","POS":"Int64"})
        knownsig = pd.concat([knownsig, knownsig_2],ignore_index=True)
        knownsig[snpid] = knownsig[snpid].astype("string")
        knownsig[snpset] = knownsig[snpset].astype("string")
        knownsig[group_key] = knownsig[group_key].astype("string")
    
    if len(knownsig)<1:
        raise ValueError("Please input a dataframe of gene list with GENE, CHR, START, END.")
    
    if group_key is not None:
        if group_key not in knownsig.columns:
            raise ValueError("Please check if group_key is in both sumstats and list of known associations.")

    ############################################################################################
    if group_key is not None:
        number_of_groups_allsig = allsig[group_key].nunique()
        number_of_groups_known = knownsig[group_key].nunique()
        log.write(" -Number of groups in sumstats:{}".format(number_of_groups_allsig), verbose=verbose)
        log.write(" -Number of groups in reference:{}".format(number_of_groups_known), verbose=verbose)

    log.write(" -Checking if variants in cis/trans regions grouped by {}...".format(group_key), verbose=verbose)
    
    ############################################################################################
    #convert to  a dict
    reference_dict = {}

    for index,row in knownsig.iterrows():
        if row[group_key] in reference_dict.keys():
            if row[snpset] in reference_dict[row[group_key]].keys():
                reference_dict[row[group_key]][row[snpset]].add(row[snpid])
            else:
                reference_dict[row[group_key]][row[snpset]] = set([row[snpid]])
        else:
            reference_dict[row[group_key]] = {row[snpset]:set([row[snpid]])}
    ############################################################################################
    #match group/trait
    try:
        no_reference_avaialble = allsig.loc[~allsig[group_key].isin(reference_dict.keys()),group_key]
        if len(no_reference_avaialble)>0:
            log.write(" -Groups not in reference: {}".format( ",".join(no_reference_avaialble)), verbose=verbose)
    except:
        pass
    ############################################################################################

    log.write(" -Checking if variants are in reference variant sets...", verbose=verbose)    
    #known_list = allsig.apply(lambda x: check_overlap(x,snpid, group_key,reference_dict), axis=1)
    new_row_list = []
    for index, row in allsig.iterrows():
        
        row = check_overlap(row, snpset, snpid, group_key,reference_dict)
        new_row_list = new_row_list+row
        known_df = pd.DataFrame(new_row_list,
                                columns=[snpid,group_key, snpset,"KNOWN_SET","OVERLAP_VARIANT","KNOWN_SET_VARIANT"])   
    
    allsig = pd.merge(allsig,known_df, on=[snpid, group_key, snpset],how="left")

    #allsig["KNOWN_SET"] = known_list.str[0]
    #allsig["OVERLAP_VARIANT"] = known_list.str[1]
    #allsig["KNOWN_SET_VARIANT"] = known_list.str[2]

    ##
    is_overlapped = ~allsig["KNOWN_SET"].isna()
    allsig["KNOWN_SET_SIZE"] = 0
    allsig.loc[is_overlapped, "KNOWN_SET_SIZE"] = allsig.loc[is_overlapped, "KNOWN_SET_VARIANT"].str.len()

    # sumstats set dic
    back_dict={}
    for i in allsig[group_key].unique():
        # for each trait in sumstats
        back_dict[i] ={}
        for j in allsig.loc[allsig[group_key]==i,snpset].unique():
            #for each locus in each trait
            back_dict[i][j] =set()
            for index, row in allsig.loc[(allsig[group_key]==i) & (allsig[snpset]==j),:].iterrows():
                #for each variant in each locus
                back_dict[i][j].add("{}".format(row["SNPID"]))

    allsig["SUMSTATS_SET_VARIANT"] = allsig.apply(lambda x: assign_set_variant(x,group_key,snpset,back_dict), axis=1)
    allsig["SUMSTATS_SET_SIZE"] = 0
    allsig["SUMSTATS_SET_SIZE"] = allsig[ "SUMSTATS_SET_VARIANT"].str.len()
    
    
    return allsig

def check_overlap(x,snpset, snpid, group_key,reference_dict):
    
    matched=[]
    if x[group_key] in reference_dict.keys():
        # if trait match
        for key, value in reference_dict[x[group_key]].items():
            # locus and snplist
            if x[snpid] in value:
                # if sumstats snp in reference snplist for locus
                # return locus and snsumstats snppid
                matched.append( (x[snpid], x[group_key],  x[snpset],  key,  x[snpid],  value))
    if len(matched)==0:
        matched = [(x[snpid], x[group_key],  x[snpset], pd.NA, pd.NA, pd.NA)]
    return matched

#def check_overlap(x,snpid, group_key,reference_dict):
#    if x[group_key] in reference_dict.keys():
#        # if trait match
#        for key, value in reference_dict[x[group_key]].items():
#            # locus and snplist
#            if x[snpid] in value:
#                # if sumstats snp in reference snplist for locus
#                # return locus and snsumstats snppid
#                return key, x[snpid], value
#    return pd.NA, pd.NA, pd.NA

def assign_set_variant(x,group_key,snpset,back_dict):
    if x[group_key] in back_dict.keys():
        # if trait match
        if x[snpset] in back_dict[x[group_key]].keys():
            #if locus match
            if len(back_dict[x[group_key]][x[snpset]]) >0:
                # return sumstats snplist for locus
                return back_dict[x[group_key]][x[snpset]]
    return pd.NA
