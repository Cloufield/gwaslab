"""
Meta-Regression for Genome-wide Association Studies (MR-MEGA)

This module implements the MR-MEGA algorithm for meta-regression of GWAS summary statistics.

Based on the methodology described in:
Mägi, R., Horikoshi, M., Sofer, T., Mahajan, A., Kitajima, H., Franceschini, N., ... & Morris, A. P. (2017).
Trans-ethnic meta-regression of genome-wide association studies accounting for ancestry increases power
for discovery and improves fine-mapping resolution. Human molecular genetics, 26(18), 3639-3650.

This is a pure Python reimplementation of the MR-MEGA algorithm, matching the original
C++ implementation exactly.

================================================================================
DETAILED WORKFLOW WITH MATHEMATICAL FORMULAS
================================================================================

The MR-MEGA algorithm performs meta-regression by:
1. Capturing genetic diversity using Multidimensional Scaling (MDS) from EAF distances
2. Using MDS coordinates as covariates in weighted least squares regression

STEP 1: Marker Selection
------------------------
Select markers that satisfy:
- Valid EAF in all cohorts: EAF_{i,j} ∈ [0, 1] for all cohorts i and markers j
- MAF threshold: min_maf < EAF_{i,j} < 1 - min_maf for all cohorts i
- Valid genomic coordinates: CHR > 0, POS > 0

From selected markers, choose one marker per 1Mb bin per chromosome (chr 1-23, bins 0-298).

STEP 2: EAF Distance Matrix Calculation
----------------------------------------
For each pair of cohorts (i, j), calculate the squared EAF distance:

    D²_{i,j} = (1/M) * Σ_{m=1}^{M} (EAF_{i,m} - EAF_{j,m})²

where:
- M = number of selected markers
- EAF_{i,m} = effect allele frequency of marker m in cohort i
- D²_{i,j} = squared distance between cohorts i and j

The distance matrix D is symmetric with D_{i,j} = D_{j,i} and D_{i,i} = 0.

STEP 3: Multidimensional Scaling (MDS)
---------------------------------------
Apply classical MDS to the distance matrix to extract principal components:

a) Double-centering transformation:
   B = -0.5 * (D - r̄ - c̄ + ḡ)
   
   where:
   - r̄_{i} = (1/n) * Σ_{j=1}^{n} D_{i,j}  (row means)
   - c̄_{j} = (1/n) * Σ_{i=1}^{n} D_{i,j}  (column means)
   - ḡ = (1/n²) * Σ_{i,j=1}^{n} D_{i,j}   (grand mean)
   - n = number of cohorts

b) Eigendecomposition:
   B = Q * Λ * Qᵀ
   
   where:
   - Q = matrix of eigenvectors (columns)
   - Λ = diagonal matrix of eigenvalues (sorted in descending order)
   - Qᵀ = transpose of Q

c) Extract MDS coordinates:
   PC_{i,k} = Q_{i,λ_k} * √(λ_k)
   
   where:
   - PC_{i,k} = k-th principal component for cohort i
   - λ_k = k-th largest eigenvalue
   - Q_{i,λ_k} = i-th element of eigenvector corresponding to λ_k

The first num_pcs components are retained as covariates.

STEP 4: Genomic Control (Optional)
-----------------------------------
For each cohort i, calculate genomic control lambda:

    λ_i = median(χ²_i) / median(χ²(1))
    
where:
- χ²_i = (β_{i,j} / SE_{i,j})² for all valid markers j in cohort i
- median(χ²(1)) = 0.4549364 (median of chi-square distribution with df=1)

If λ_i > 1, adjust standard errors:
    SE_{i,j}' = SE_{i,j} * √(λ_i)

STEP 5: Weighted Least Squares Regression
-----------------------------------------
For each variant, perform weighted least squares regression:

    Y = X * β + ε

where:
- Y = vector of effect sizes (β values) from valid cohorts, shape (n_valid,)
- X = design matrix [1, PC₁, PC₂, ..., PC_{num_pcs}], shape (n_valid, num_pcs+1)
- β = regression coefficients [β₀, β₁, β₂, ..., β_{num_pcs}], shape (num_pcs+1,)
- ε = error term
- n_valid = number of cohorts with valid data for this variant

The design matrix X is:
    X = [1  PC₁₁  PC₁₂  ...  PC₁_{num_pcs}]
        [1  PC₂₁  PC₂₂  ...  PC₂_{num_pcs}]
        [...                              ]
        [1  PC_{n_valid,1}  ...  PC_{n_valid,num_pcs}]

Weighted least squares solution:
    β̂ = (Xᵀ W X)⁻¹ Xᵀ W Y

where:
- W = diagonal weight matrix with W_{i,i} = 1 / SE²_i
- SE_i = standard error for cohort i (after genomic control if applied)

Variance-covariance matrix:
    Var(β̂) = MSE * (Xᵀ W X)⁻¹

where:
    MSE = RSS / (n_valid - num_pcs - 1)
    RSS = Σ_{i=1}^{n_valid} W_{i,i} * (Y_i - Ŷ_i)²  (residual sum of squares)
    Ŷ = X * β̂  (predicted values)

Standard errors:
    SE(β̂_k) = √(Var(β̂)_{k,k})

STEP 6: Statistical Tests
--------------------------
a) Association test (TSS0 - RSS):
    χ²_assoc = |TSS0 - RSS|
    df_assoc = num_pcs + 1
    
    where:
    - TSS0 = Σ_{i=1}^{n_valid} W_{i,i} * Y²_i  (total sum of squares, centered at 0)
    - RSS = Σ_{i=1}^{n_valid} W_{i,i} * (Y_i - Ŷ_i)²  (residual sum of squares)
    
    P-value: P_assoc = 1 - F_χ²(χ²_assoc, df_assoc)

b) Ancestry heterogeneity test (TSS - RSS):
    χ²_ancestry = |TSS - RSS|
    df_ancestry = num_pcs
    
    where:
    - TSS = Σ_{i=1}^{n_valid} W_{i,i} * (Y_i - Ȳ)²  (total sum of squares, centered at mean)
    - Ȳ = (Σ_{i=1}^{n_valid} W_{i,i} * Y_i) / (Σ_{i=1}^{n_valid} W_{i,i})  (weighted mean)
    
    P-value: P_ancestry = 1 - F_χ²(χ²_ancestry, df_ancestry)

c) Residual heterogeneity test (RSS):
    χ²_residual = |RSS|
    df_residual = n_valid - num_pcs - 1
    
    P-value: P_residual = 1 - F_χ²(χ²_residual, df_residual)

d) Bayes factor:
    ln(BF) = (TSS0_RSS - (num_pcs + 1) * ln(n_valid)) / 2
    
    where TSS0_RSS = |TSS0 - RSS|

STEP 7: Second Genomic Control Correction (GCO, Optional)
---------------------------------------------------------
If use_gco is enabled, apply a second genomic control correction:

a) Calculate GCO lambda:
    λ_gco = median(χ²_assoc) / median(χ²(df_assoc))
    
    where:
    - median(χ²_assoc) = median of all valid association chi-squares
    - median(χ²(df_assoc)) = median of chi-square distribution with df = df_assoc

b) Correct association chi-squares:
    χ²_assoc' = χ²_assoc / λ_gco

c) Recalculate P-values:
    P_assoc' = 1 - F_χ²(χ²_assoc', df_assoc)

OUTPUT
------
For each variant, the algorithm outputs:
- β₀, SE(β₀): Intercept (main association effect)
- β₁, SE(β₁), ..., β_{num_pcs}, SE(β_{num_pcs}): PC coefficients
- χ²_assoc, df_assoc, P_assoc: Association test statistics
- χ²_ancestry, df_ancestry, P_ancestry: Ancestry heterogeneity test statistics
- χ²_residual, df_residual, P_residual: Residual heterogeneity test statistics
- ln(BF): Natural logarithm of Bayes factor
- EAF: Weighted average EAF across cohorts
- Nsample: Total sample size
- Ncohort: Number of cohorts with valid data
- Effects: Direction of effects across cohorts (+/-/0/?)

================================================================================
"""

from typing import TYPE_CHECKING, Optional
import numpy as np
import pandas as pd
import os
import subprocess
import shutil
import tempfile
from scipy.stats import norm
from scipy.stats.distributions import chi2
from scipy.linalg import inv, pinv, eigh, solve
from gwaslab.info.g_Log import Log
from gwaslab.g_Sumstats import Sumstats

if TYPE_CHECKING:
    from gwaslab.g_SumstatsMulti import SumstatsMulti


def _parse_mrmega_output(
    output_file: str,
    log: Optional[Log] = None,
    verbose: bool = True
) -> Optional[Sumstats]:
    """
    Parse MR-MEGA output file (.result file).
    
    MR-MEGA output columns:
    - MarkerName, Chromosome, Position, EA, NEA, EAF, Nsample, Ncohort, Effects
    - beta_0, se_0, beta_1, se_1, ... (PC effects)
    - chisq_association, ndf_association, P-value_association
    - chisq_ancestry_het, ndf_ancestry_het, P-value_ancestry_het
    - chisq_residual_het, ndf_residual_het, P-value_residual_het
    - lnBF, Comments
    
    Parameters
    ----------
    output_file : str
        Path to MR-MEGA .result file.
    log : Log, optional
        Log object for logging messages.
    verbose : bool, default True
        Whether to print progress messages.
    
    Returns
    -------
    Sumstats or None
        Parsed MR-MEGA results as Sumstats object.
    """
    if log is None:
        log = Log()
    
    try:
        # Read MR-MEGA result file (tab-delimited)
        result_df = pd.read_csv(output_file, sep="\t", comment="#")
        
        log.write(f"MR-MEGA output columns: {list(result_df.columns)}", verbose=verbose)
        log.write(f"MR-MEGA output shape: {result_df.shape}", verbose=verbose)
        
        # Map MR-MEGA columns to GWASLab format
        column_mapping = {
            "MarkerName": "SNPID",
            "Chromosome": "CHR",
            "Position": "POS",
            "EA": "EA",
            "NEA": "NEA",
            "EAF": "EAF",
            "Nsample": "N",
        }
        
        # Rename columns
        for old_col, new_col in column_mapping.items():
            if old_col in result_df.columns and new_col not in result_df.columns:
                result_df = result_df.rename(columns={old_col: new_col})
        
        # Extract main effect (beta_0 is the intercept in MR-MEGA's regression model)
        # MR-MEGA model: Y = beta_0 + beta_1*PC1 + beta_2*PC2 + ... + error
        # beta_0 is the intercept (main association effect)
        # beta_1, beta_2, ... are the PC coefficients
        # Handle "NA" strings as well as actual NaN values
        if "beta_0" in result_df.columns:
            # Convert "NA" strings to NaN
            result_df["beta_0"] = result_df["beta_0"].replace("NA", np.nan)
            result_df["BETA"] = pd.to_numeric(result_df["beta_0"], errors='coerce')
        if "se_0" in result_df.columns:
            result_df["se_0"] = result_df["se_0"].replace("NA", np.nan)
            result_df["SE"] = pd.to_numeric(result_df["se_0"], errors='coerce')
        
        # Use P-value_association as the main p-value
        if "P-value_association" in result_df.columns:
            result_df["P-value_association"] = result_df["P-value_association"].replace("NA", np.nan)
            result_df["P"] = pd.to_numeric(result_df["P-value_association"], errors='coerce')
        
        # Calculate Z-score from BETA and SE if available
        if "BETA" in result_df.columns and "SE" in result_df.columns:
            result_df["Z"] = result_df["BETA"] / result_df["SE"]
        
        # Map heterogeneity statistics (handle "NA" strings)
        if "chisq_residual_het" in result_df.columns:
            result_df["chisq_residual_het"] = result_df["chisq_residual_het"].replace("NA", np.nan)
            result_df["Q"] = pd.to_numeric(result_df["chisq_residual_het"], errors='coerce')
        if "ndf_residual_het" in result_df.columns:
            result_df["ndf_residual_het"] = result_df["ndf_residual_het"].replace("NA", np.nan)
            result_df["DOF"] = pd.to_numeric(result_df["ndf_residual_het"], errors='coerce').astype("Int64")
        if "P-value_residual_het" in result_df.columns:
            result_df["P-value_residual_het"] = result_df["P-value_residual_het"].replace("NA", np.nan)
            result_df["P_HET"] = pd.to_numeric(result_df["P-value_residual_het"], errors='coerce')
        
        # Calculate I2 if Q and DOF are available
        if "Q" in result_df.columns and "DOF" in result_df.columns:
            # I2 = max(0, (Q - df) / Q)
            result_df["I2"] = np.where(
                result_df["Q"] > 0,
                np.maximum(0.0, (result_df["Q"] - result_df["DOF"]) / result_df["Q"]),
                0.0
            )
        
        # Ensure required columns exist
        required_cols = ["SNPID", "CHR", "POS", "EA", "NEA"]
        if not all(col in result_df.columns for col in required_cols):
            log.write(
                f"MR-MEGA output missing required columns. Found: {list(result_df.columns)}",
                verbose=verbose
            )
            return None
        
        # Create Sumstats object
        # Keep MR-MEGA specific columns as "other" columns
        standard_cols = required_cols + ["BETA", "SE", "P", "Z", "EAF", "N", "Q", "DOF", "P_HET", "I2"]
        other_cols = [col for col in result_df.columns if col not in standard_cols]
        
        result_sumstats = Sumstats(result_df, fmt="gwaslab", other=other_cols)
        return result_sumstats
        
    except Exception as e:
        log.write(f"Error parsing MR-MEGA output: {str(e)}", verbose=verbose)
        import traceback
        log.write(traceback.format_exc(), verbose=verbose)
        return None


def _check_mrmega_available(mrmega_path: Optional[str] = None) -> Optional[str]:
    """
    Check if MR-MEGA is available in PATH or at specified path.
    
    Parameters
    ----------
    mrmega_path : str, optional
        Path to MR-MEGA executable. If None, searches PATH.
    
    Returns
    -------
    str or None
        Path to MR-MEGA executable if found, None otherwise.
    """
    if mrmega_path is not None:
        if os.path.isfile(mrmega_path) and os.access(mrmega_path, os.X_OK):
            return mrmega_path
        return None
    
    # Search in PATH
    mrmega_path = shutil.which("MR-MEGA") or shutil.which("mrmega") or shutil.which("MRMEGA")
    return mrmega_path


def _calculate_eaf_distance_matrix(
    eaf_matrix: np.ndarray
) -> np.ndarray:
    """
    Calculate distance matrix from EAF (Effect Allele Frequency) values.
    
    Distance between cohorts i and j is: mean((EAF_k_i - EAF_k_j)^2) over all markers k.
    This matches MR-MEGA's distance() method in structures.cpp.
    
    Parameters
    ----------
    eaf_matrix : np.ndarray
        EAF matrix of shape (n_markers, n_cohorts)
    
    Returns
    -------
    np.ndarray
        Distance matrix of shape (n_cohorts, n_cohorts)
    """
    n_markers, n_cohorts = eaf_matrix.shape
    # eaf_matrix should already be float64 from caller, but ensure for safety
    if eaf_matrix.dtype != np.float64:
        eaf_matrix = eaf_matrix.astype(np.float64, copy=False)
    
    # Vectorized distance calculation: use broadcasting to compute all pairwise differences
    # eaf_matrix[:, :, None] - eaf_matrix[:, None, :] gives (n_markers, n_cohorts, n_cohorts)
    # Then sum over markers and divide by n_markers
    diff_sq = np.sum((eaf_matrix[:, :, None] - eaf_matrix[:, None, :]) ** 2, axis=0, dtype=np.float64)
    dist_matrix = diff_sq / n_markers
    
    return dist_matrix


def _calculate_mds_from_distance(
    dist_matrix: np.ndarray,
    n_components: int
) -> np.ndarray:
    """
    Calculate MDS (Multidimensional Scaling) from distance matrix.
    
    This implements MR-MEGA's calculateMDS function using double-centering
    and eigendecomposition.
    
    Parameters
    ----------
    dist_matrix : np.ndarray
        Distance matrix of shape (n_cohorts, n_cohorts)
    n_components : int
        Number of MDS components to extract
    
    Returns
    -------
    np.ndarray
        MDS coordinates of shape (n_cohorts, n_components)
    """
    n_cohorts = dist_matrix.shape[0]
    
    # dist_matrix should already be float64 from caller, but ensure for safety
    if dist_matrix.dtype != np.float64:
        dist_matrix = dist_matrix.astype(np.float64, copy=False)
    
    # Double-centering: B = -0.5 * (D - row_means - col_means + grand_mean)
    # This is equivalent to: B[i,j] = -0.5 * (D[i,j] - mean_i - mean_j + mean_all)
    row_means = dist_matrix.mean(axis=1, keepdims=True, dtype=np.float64)
    col_means = dist_matrix.mean(axis=0, keepdims=True, dtype=np.float64)
    grand_mean = dist_matrix.mean(dtype=np.float64)
    
    # Double-centered matrix (result is float64)
    B = -0.5 * (dist_matrix - row_means - col_means + grand_mean)
    
    # Eigendecomposition (for symmetric matrix, equivalent to SVD)
    # MR-MEGA uses SVD, but for symmetric matrices, eigendecomposition gives same result
    # eigh returns float64 by default
    eigenvalues, eigenvectors = eigh(B)
    
    # Sort by eigenvalue (descending) and get indices
    idx = np.argsort(eigenvalues)[::-1]
    eigenvalues_sorted = eigenvalues[idx]
    eigenvectors_sorted = eigenvectors[:, idx]
    
    # Take top n_components
    n_components = min(n_components, n_cohorts)
    
    # Create list of eigenvalue indices (matching MR-MEGA's elist)
    elist = idx[:n_components].tolist()
    
    # Calculate sqrt of eigenvalues (only non-negative) - vectorized
    # Result is float64 from np.where and np.sqrt
    sqrt_eigenvalues = np.where(eigenvalues >= 0, np.sqrt(eigenvalues), 0.0)
    
    # Vectorized MDS coordinate calculation
    mds_coords = np.zeros((n_cohorts, n_components), dtype=np.float64)
    for c2 in range(n_components):  # PC index
        i2 = elist[c2]  # original eigenvalue index (before sorting)
        # Vectorize across cohorts: multiply eigenvector column by sqrt_eigenvalue
        mds_coords[:, c2] = eigenvectors_sorted[:, c2] * sqrt_eigenvalues[i2]
    
    return mds_coords


def _weighted_least_squares(
    Y: np.ndarray,
    X: np.ndarray,
    W: np.ndarray
) -> tuple:
    """
    Perform weighted least squares regression.
    
    This implements MR-MEGA's lr_w function from regression.cpp.
    
    Model: Y = X @ beta + error
    Weights: W (inverse variance)
    
    Parameters
    ----------
    Y : np.ndarray
        Dependent variable (beta values), shape (n,)
    X : np.ndarray
        Design matrix, shape (n, p) where p is number of predictors
    W : np.ndarray
        Weights (inverse variance), shape (n,)
    
    Returns
    -------
    tuple
        (beta_hat, se_beta, TSS, TSS0, RSS, var_cov)
        - beta_hat: coefficient estimates, shape (p,)
        - se_beta: standard errors, shape (p,)
        - TSS: total sum of squares (weighted, centered at mean)
        - TSS0: total sum of squares (weighted, centered at 0)
        - RSS: residual sum of squares (weighted)
        - var_cov: variance-covariance matrix, shape (p, p)
    """
    n, p = X.shape
    df = n - p
    
    if df < 1:
        return None, None, None, None, None, None
    
    # Inputs should already be float64, but ensure for safety (only if needed)
    if X.dtype != np.float64:
        X = X.astype(np.float64, copy=False)
    if Y.dtype != np.float64:
        Y = Y.astype(np.float64, copy=False)
    if W.dtype != np.float64:
        W = W.astype(np.float64, copy=False)
    
    # Use einsum for efficient weighted matrix operations (avoids intermediate arrays)
    XtWX = np.einsum('ki,k,kj->ij', X, W, X, optimize=True)
    
    # Reuse W * Y computation (needed for both XtWy and Y_bar/TSS0)
    WY = W * Y
    XtWy = np.einsum('ki,k->i', X, WY, optimize=True)
    
    # Use solve instead of inversion (faster and more numerically stable)
    try:
        beta_hat = np.linalg.solve(XtWX, XtWy)
        # Only compute inverse if needed for variance-covariance matrix
        try:
            XtWX_inv = inv(XtWX)
        except np.linalg.LinAlgError:
            XtWX_inv = pinv(XtWX)
    except np.linalg.LinAlgError:
        # Fall back to pseudo-inverse if singular
        XtWX_inv = pinv(XtWX)
        beta_hat = XtWX_inv @ XtWy
    
    # Calculate predicted values and residuals
    Y_pred = X @ beta_hat
    residuals = Y - Y_pred
    
    # Calculate weighted mean of Y (reuse WY)
    W_sum = np.sum(W, dtype=np.float64)
    Y_bar = np.sum(WY, dtype=np.float64) / W_sum
    
    # Optimized sums of squares (reuse WY and compute Y-Y_bar once)
    Y_centered = Y - Y_bar
    TSS = np.sum(W * Y_centered ** 2, dtype=np.float64)
    TSS0 = np.sum(WY * Y, dtype=np.float64)  # Reuse WY: W * Y^2 = (W*Y) * Y
    RSS = np.sum(W * residuals ** 2, dtype=np.float64)
    
    # Mean squared error (result is float64)
    MSE = RSS / df if df > 0 else 0.0
    
    # Variance-covariance matrix: MSE * (X' W X)^(-1) (result is float64)
    var_cov = MSE * XtWX_inv
    
    # Standard errors (result is float64)
    se_beta = np.sqrt(np.diag(var_cov))
    
    return beta_hat, se_beta, TSS, TSS0, RSS, var_cov


def _weighted_least_squares_batch(
    Y: np.ndarray,
    X: np.ndarray,
    W: np.ndarray
) -> tuple:
    """
    Vectorized batch version of weighted least squares regression.
    
    Processes multiple variants with the same valid study pattern simultaneously.
    Uses SciPy's batched linear algebra operations for improved performance.
    See: https://docs.scipy.org/doc/scipy-1.17.0/tutorial/linalg_batch.html
    
    Model: Y = X @ beta + error (for each variant)
    Weights: W (inverse variance)
    
    Parameters
    ----------
    Y : np.ndarray
        Dependent variable (beta values), shape (n_variants, n_valid)
    X : np.ndarray
        Design matrix, shape (n_valid, p) where p is number of predictors
    W : np.ndarray
        Weights (inverse variance), shape (n_variants, n_valid)
    
    Returns
    -------
    tuple
        (beta_hat, se_beta, TSS, TSS0, RSS)
        - beta_hat: coefficient estimates, shape (n_variants, p)
        - se_beta: standard errors, shape (n_variants, p)
        - TSS: total sum of squares (weighted, centered at mean), shape (n_variants,)
        - TSS0: total sum of squares (weighted, centered at 0), shape (n_variants,)
        - RSS: residual sum of squares (weighted), shape (n_variants,)
    
    Notes
    -----
    This function uses SciPy's batched linear algebra operations (solve, inv) which
    process all variants simultaneously for better performance. If batched operations
    fail (e.g., due to singular matrices), it falls back to per-variant processing.
    """
    n_variants, n_valid = Y.shape
    p = X.shape[1]
    df = n_valid - p
    
    if df < 1:
        return None, None, None, None, None
    
    # Ensure float64
    if X.dtype != np.float64:
        X = X.astype(np.float64, copy=False)
    if Y.dtype != np.float64:
        Y = Y.astype(np.float64, copy=False)
    if W.dtype != np.float64:
        W = W.astype(np.float64, copy=False)
    
    # Pre-compute XtWX once (same for all variants since X is the same)
    # XtWX = X^T @ diag(W_mean) @ X, but we need per-variant weights
    # Actually, we need to compute XtWX for each variant since W varies
    # Use einsum for batch computation: 'vki,k,vk,vkj->vij'
    # where v = variant index, k = study index, i/j = predictor index
    XtWX = np.einsum('ki,vk,kj->vij', X, W, X, optimize=True)  # (n_variants, p, p)
    
    # Compute WY for all variants
    WY = W * Y  # (n_variants, n_valid)
    
    # Compute XtWy for all variants: X^T @ (W * Y)
    XtWy = np.einsum('ki,vk->vi', X, WY, optimize=True)  # (n_variants, p)
    
    # Use batched linear algebra operations from SciPy
    # According to https://docs.scipy.org/doc/scipy-1.17.0/tutorial/linalg_batch.html
    # SciPy's solve() and inv() support batched operations when input arrays have
    # batch dimensions (first dimensions) followed by core dimensions (last 2 dims)
    beta_hat = np.zeros((n_variants, p), dtype=np.float64)
    XtWX_inv = np.zeros((n_variants, p, p), dtype=np.float64)
    
    # Reshape XtWy for batched solve: (n_variants, p) -> (n_variants, p, 1)
    # This matches the batched solve API: solve(A, b) where A is (batch, n, n) and b is (batch, n, 1)
    XtWy_batched = XtWy[:, :, np.newaxis]  # (n_variants, p, 1)
    
    try:
        # Attempt batched solve: solve all linear systems at once
        # XtWX shape: (n_variants, p, p) - batch shape (n_variants,), core shape (p, p)
        # XtWy_batched shape: (n_variants, p, 1) - batch shape (n_variants,), core shape (p, 1)
        # Result shape: (n_variants, p, 1)
        beta_hat_batched = solve(XtWX, XtWy_batched)
        beta_hat = beta_hat_batched[:, :, 0]  # Extract from (n_variants, p, 1) to (n_variants, p)
        
        # Attempt batched inverse: compute all matrix inverses at once
        # XtWX shape: (n_variants, p, p) - batch shape (n_variants,), core shape (p, p)
        # Result shape: (n_variants, p, p)
        XtWX_inv = inv(XtWX)
        
    except (np.linalg.LinAlgError, ValueError) as e:
        # Fall back to per-variant processing if batched operations fail
        # This can happen if some matrices in the batch are singular
        # Process each variant individually with proper error handling
        for v in range(n_variants):
            try:
                beta_hat[v] = solve(XtWX[v], XtWy[v])
                try:
                    XtWX_inv[v] = inv(XtWX[v])
                except np.linalg.LinAlgError:
                    XtWX_inv[v] = pinv(XtWX[v])
            except np.linalg.LinAlgError:
                # Use pseudo-inverse for singular matrices
                XtWX_inv_v = pinv(XtWX[v])
                XtWX_inv[v] = XtWX_inv_v
                beta_hat[v] = XtWX_inv_v @ XtWy[v]
    
    # Calculate predicted values and residuals for all variants
    Y_pred = (X @ beta_hat.T).T  # (n_variants, n_valid)
    residuals = Y - Y_pred  # (n_variants, n_valid)
    
    # Calculate weighted mean of Y for each variant
    W_sum = np.sum(W, axis=1, dtype=np.float64, keepdims=True)  # (n_variants, 1)
    Y_bar = np.sum(WY, axis=1, dtype=np.float64, keepdims=True) / W_sum  # (n_variants, 1)
    
    # Calculate sums of squares for all variants
    Y_centered = Y - Y_bar  # (n_variants, n_valid)
    TSS = np.sum(W * Y_centered ** 2, axis=1, dtype=np.float64)  # (n_variants,)
    TSS0 = np.sum(WY * Y, axis=1, dtype=np.float64)  # (n_variants,)
    RSS = np.sum(W * residuals ** 2, axis=1, dtype=np.float64)  # (n_variants,)
    
    # Calculate MSE and standard errors
    MSE = RSS / df if df > 0 else np.zeros(n_variants, dtype=np.float64)
    
    # Variance-covariance matrices: MSE * (X' W X)^(-1)
    # Extract diagonal (standard errors) efficiently
    var_diag = np.einsum('v,vii->vi', MSE, XtWX_inv)  # (n_variants, p)
    se_beta = np.sqrt(var_diag)  # (n_variants, p)
    
    return beta_hat, se_beta, TSS, TSS0, RSS


def meta_regress_mrmega_python(
    sumstats_multi: 'SumstatsMulti',
    num_pcs: int = 2,
    use_genomic_control: bool = True,
    use_gco: bool = True,
    min_maf: float = 0.01,
    log: Optional[Log] = None,
    verbose: bool = True,
    selected_marker_snpids: Optional[list] = None
) -> Optional[Sumstats]:
    """
    Perform meta-regression using MR-MEGA algorithm implemented in pure Python.
    
    This reimplements the MR-MEGA algorithm without requiring the external executable.
    The algorithm:
    
    1. Selects markers with MAF > min_maf in all cohorts
    2. Calculates EAF correlation distance matrix between cohorts
    3. Performs MDS to extract principal components
    4. For each marker, performs weighted least squares regression:
       Y = beta_0 + beta_1*PC1 + beta_2*PC2 + ... + error
    5. Calculates association, ancestry heterogeneity, and residual heterogeneity statistics
    
    Parameters
    ----------
    sumstats_multi : SumstatsMulti
        SumstatsMulti object containing summary statistics from multiple studies.
    num_pcs : int, default 2
        Number of principal components to use. Must satisfy: num_pcs < nstudy - 2
    use_genomic_control : bool, default True
        Whether to apply genomic control correction to standard errors (--gc flag).
    use_gco : bool, default True
        Whether to apply second genomic control correction to output chi-squares (--gco flag).
    min_maf : float, default 0.01
        Minimum minor allele frequency for marker selection (1%).
    log : Log, optional
        Log object for logging messages. If None, creates a new Log.
    verbose : bool, default True
        Whether to print progress messages.
    
    Returns
    -------
    Sumstats or None
        Sumstats object containing MR-MEGA results.
    """
    if log is None:
        log = Log()
    
    nstudy = sumstats_multi.meta["gwaslab"]["number_of_studies"]
    
    # Check constraint: num_pcs < nstudy - 2
    max_allowed_pcs = max(1, nstudy - 3)
    if num_pcs > max_allowed_pcs:
        log.write(
            f"Warning: Requested {num_pcs} PCs but max allowed is {max_allowed_pcs} "
            f"(nstudy={nstudy}, constraint: npcs < nstudy-2). Using {max_allowed_pcs} PCs instead.",
            verbose=verbose
        )
        num_pcs = max_allowed_pcs
    
    log.write(f"Running MR-MEGA algorithm in Python (nstudy={nstudy}, npcs={num_pcs})", verbose=verbose)
    log.write(
        "Citation: Mägi, R., Horikoshi, M., Sofer, T., Mahajan, A., Kitajima, H., "
        "Franceschini, N., ... & Morris, A. P. (2017). Trans-ethnic meta-regression of "
        "genome-wide association studies accounting for ancestry increases power for discovery "
        "and improves fine-mapping resolution. Human molecular genetics, 26(18), 3639-3650.",
        verbose=verbose
    )
    
    # Get data
    if sumstats_multi.engine == "polars":
        data = sumstats_multi.data.to_pandas()
    else:
        data = sumstats_multi.data.copy()
    
    data = data.sort_values(by=["CHR", "POS"]).reset_index(drop=True)
    
    # Extract columns
    beta_cols = [f"BETA_{i+1}" for i in range(nstudy)]
    se_cols = [f"SE_{i+1}" for i in range(nstudy)]
    n_cols = [f"N_{i+1}" for i in range(nstudy)]
    eaf_cols = [f"EAF_{i+1}" for i in range(nstudy)]
    
    # Step 1: Select markers for MDS calculation
    # Markers must have:
    # - Valid EAF in all cohorts
    # - MAF > min_maf in all cohorts (0.01 < EAF < 0.99)
    # - Valid CHR and POS
    log.write("Selecting markers for MDS calculation...", verbose=verbose)
    
    use_pre_selected = False
    if selected_marker_snpids is not None and len(selected_marker_snpids) > 0:
        selected_mask = data["SNPID"].isin(selected_marker_snpids)
        selected_markers = data.loc[selected_mask].copy()
        selected_indices = selected_markers.index
        
        if len(selected_markers) == 0:
            log.write("Warning: None of the pre-selected markers found in data. Using default selection.", verbose=verbose)
            use_pre_selected = False
        else:
            use_pre_selected = True
    
    if not use_pre_selected:
        valid_eaf_mask = data[eaf_cols].notna().all(axis=1)
        valid_chr_pos = data["CHR"].notna() & data["POS"].notna() & (data["CHR"] > 0) & (data["POS"] > 0)
        
        maf_mask = (
            (data[eaf_cols] > min_maf).all(axis=1) &
            (data[eaf_cols] < (1 - min_maf)).all(axis=1)
        )
        
        good_markers_mask = valid_eaf_mask & valid_chr_pos & maf_mask
        
        if good_markers_mask.sum() > 0:
            good_markers = data.loc[good_markers_mask].copy()
            good_markers["BIN"] = good_markers["POS"] // 1000000
            
            good_markers = good_markers.reset_index(drop=True)
            
            matrix = {}
            markerNum = 1
            
            for idx, row in good_markers.iterrows():
                chr = int(row["CHR"])
                bin_num = int(row["BIN"])
                matrix[(chr, bin_num)] = markerNum - 1
                markerNum += 1
            
            selected_indices_list = []
            for chr in range(1, 24):
                for bin_num in range(0, 299):
                    if (chr, bin_num) in matrix:
                        stored_value = matrix[(chr, bin_num)]
                        if stored_value > 0:
                            selected_indices_list.append(stored_value)
            
            if len(selected_indices_list) > 0:
                selected_indices = good_markers.index[selected_indices_list]
                selected_markers = good_markers.loc[selected_indices].copy().reset_index(drop=True)
            else:
                selected_markers = pd.DataFrame(columns=good_markers.columns)
                selected_indices = pd.Index([])
            
            log.write(
                f"Selected {len(selected_markers)} markers for MDS calculation "
                f"(from {good_markers_mask.sum()} markers with MAF >= {min_maf})",
                verbose=verbose
            )
        else:
            log.write(
                f"Warning: No markers found with MAF > {min_maf} in all cohorts. "
                f"Using all available markers.",
                verbose=verbose
            )
            selected_markers = data.loc[valid_eaf_mask].copy()
            selected_indices = selected_markers.index
    
    # Step 2: Calculate distance matrix from EAF
    log.write("Calculating EAF distance matrix...", verbose=verbose)
    eaf_matrix = selected_markers[eaf_cols].values  # shape: (n_markers, n_cohorts)
    dist_matrix = _calculate_eaf_distance_matrix(eaf_matrix)
    
    # Step 3: Calculate MDS
    log.write(f"Calculating MDS with {num_pcs} components...", verbose=verbose)
    pcs = _calculate_mds_from_distance(dist_matrix, num_pcs)  # shape: (n_cohorts, num_pcs)
    
    log.write("Principal components calculated:", verbose=verbose)
    for i in range(nstudy):
        pc_str = " ".join([f"{pcs[i, j]:.6f}" for j in range(num_pcs)])
        log.write(f"  Cohort {i+1}: {pc_str}", verbose=verbose)
    
    # Step 4: Calculate genomic control lambda for each cohort (if enabled)
    # Pre-extract data arrays for GC calculation (before main loop)
    lambda_cohorts = np.ones(nstudy, dtype=np.float64)
    if use_genomic_control:
        log.write("Calculating genomic control lambda for each cohort...", verbose=verbose)
        # Use pre-extracted arrays for faster computation
        beta_gc = data[beta_cols].values.astype(np.float64, copy=False)
        se_gc = data[se_cols].values.astype(np.float64, copy=False)
        for i in range(nstudy):
            valid = ~(np.isnan(beta_gc[:, i]) | np.isnan(se_gc[:, i]) | (se_gc[:, i] <= 0))
            if valid.sum() > 0:
                chi2_values = (beta_gc[valid, i] / se_gc[valid, i]) ** 2
                median_chi2 = np.median(chi2_values)
                # Lambda = median_chi2 / median_of_chi2(df=1)
                lambda_cohorts[i] = median_chi2 / 0.4549364  # median of chi2(1) distribution
                if lambda_cohorts[i] > 1:
                    log.write(f"  Cohort {i+1}: lambda = {lambda_cohorts[i]:.4f}", verbose=verbose)
    
    # Step 5: Perform regression for each marker
    log.write("Performing meta-regression for each marker...", verbose=verbose)
    
    n_variants = len(data)
    # Initialize results arrays with float64 dtype to avoid conversions
    results = {
        "beta_0": np.full(n_variants, np.nan, dtype=np.float64),
        "se_0": np.full(n_variants, np.nan, dtype=np.float64),
    }
    # Add PC coefficients
    for j in range(num_pcs):
        results[f"beta_{j+1}"] = np.full(n_variants, np.nan, dtype=np.float64)
        results[f"se_{j+1}"] = np.full(n_variants, np.nan, dtype=np.float64)
    
    results.update({
        "chisq_association": np.full(n_variants, np.nan, dtype=np.float64),
        "ndf_association": np.full(n_variants, np.nan, dtype=np.float64),
        "P-value_association": np.full(n_variants, np.nan, dtype=np.float64),
        "chisq_ancestry_het": np.full(n_variants, np.nan, dtype=np.float64),
        "ndf_ancestry_het": np.full(n_variants, np.nan, dtype=np.float64),
        "P-value_ancestry_het": np.full(n_variants, np.nan, dtype=np.float64),
        "chisq_residual_het": np.full(n_variants, np.nan, dtype=np.float64),
        "ndf_residual_het": np.full(n_variants, np.nan, dtype=np.float64),
        "P-value_residual_het": np.full(n_variants, np.nan, dtype=np.float64),
        "lnBF": np.full(n_variants, np.nan, dtype=np.float64),
        "Comments": [""] * n_variants,
        "Ncohort": np.zeros(n_variants, dtype=int),
        "Nsample": np.zeros(n_variants, dtype=int),
        "EAF": np.full(n_variants, np.nan, dtype=np.float64),
        "Effects": [""] * n_variants,
    })
    
    # Store association chi-squares for second GC correction (gco)
    # MR-MEGA has --gco option for second genomic control on output
    # Pre-allocate as numpy array for better performance
    good_chisq_assoc = np.full(n_variants, np.nan, dtype=np.float64)
    good_chisq_idx = 0
    
    # Pre-extract data arrays for faster access (avoid repeated iloc calls)
    # Ensure float64 at extraction to avoid repeated conversions
    beta_data = data[beta_cols].values.astype(np.float64, copy=False)  # (n_variants, nstudy)
    se_data = data[se_cols].values.astype(np.float64, copy=False)  # (n_variants, nstudy)
    n_data = data[n_cols].values.astype(np.float64, copy=False) if n_cols[0] in data.columns else None
    eaf_data = data[eaf_cols].values.astype(np.float64, copy=False) if eaf_cols[0] in data.columns else None
    
    snpids = data["SNPID"].values if "SNPID" in data.columns else np.array([f"variant_{i}" for i in range(n_variants)])
    
    # Pre-compute valid masks for all variants (vectorized)
    # valid_masks[i, j] = True if variant i has valid data for cohort j
    valid_masks = ~(np.isnan(beta_data) | np.isnan(se_data) | (se_data <= 0))  # (n_variants, nstudy)
    
    # Pre-compute count of valid studies per variant (vectorized)
    n_valid_per_variant = np.count_nonzero(valid_masks, axis=1)  # (n_variants,)
    
    # Pre-apply genomic control to SE data if enabled (vectorized)
    if use_genomic_control:
        # Use broadcasting: lambda_cohorts[None, :] broadcasts to (n_variants, nstudy)
        gc_mask_matrix = lambda_cohorts[None, :] > 1  # (1, nstudy) -> (n_variants, nstudy)
        # Apply GC correction: SE *= sqrt(lambda) where lambda > 1
        # Broadcasting: sqrt(lambda_cohorts[None, :]) -> (1, nstudy) -> (n_variants, nstudy)
        se_data = np.where(gc_mask_matrix, se_data * np.sqrt(lambda_cohorts[None, :]), se_data)
    
    # Cache constants and function references
    num_pcs_plus_1 = num_pcs + 1
    chi2_sf = chi2.sf  # Cache function reference to reduce lookup overhead
    
    # Cache result arrays to reduce dictionary lookups
    result_beta_0 = results["beta_0"]
    result_se_0 = results["se_0"]
    result_beta_pc = [results[f"beta_{j+1}"] for j in range(num_pcs)]
    result_se_pc = [results[f"se_{j+1}"] for j in range(num_pcs)]
    result_chisq_assoc = results["chisq_association"]
    result_ndf_assoc = results["ndf_association"]
    result_p_assoc = results["P-value_association"]
    result_chisq_ancestry = results["chisq_ancestry_het"]
    result_ndf_ancestry = results["ndf_ancestry_het"]
    result_p_ancestry = results["P-value_ancestry_het"]
    result_chisq_residual = results["chisq_residual_het"]
    result_ndf_residual = results["ndf_residual_het"]
    result_p_residual = results["P-value_residual_het"]
    result_lnBF = results["lnBF"]
    
    # Group variants by valid study pattern for batch processing
    # Convert valid_masks to tuples for hashing
    # Create pattern groups: variants with same valid study pattern can be processed together
    # Use a dictionary to group variants by their valid study pattern
    pattern_groups = {}
    for idx in range(n_variants):
        n_valid = n_valid_per_variant[idx]
        if n_valid == 0:
            results["Comments"][idx] = "NoValidStudies"
            continue
        if n_valid - 2 <= num_pcs:
            results["Comments"][idx] = "SmallCohortCount"
            continue
        
        # Create a tuple key from the valid mask (hashable)
        valid_mask = valid_masks[idx, :]
        pattern_key = tuple(valid_mask)
        
        if pattern_key not in pattern_groups:
            pattern_groups[pattern_key] = []
        pattern_groups[pattern_key].append(idx)
    
    # Process each pattern group in batches
    for pattern_key, variant_indices in pattern_groups.items():
        if len(variant_indices) == 0:
            continue
        
        # Convert variant_indices to numpy array for efficient indexing
        variant_indices = np.array(variant_indices, dtype=int)
        
        # Get valid studies for this pattern (same for all variants in group)
        valid_mask = np.array(pattern_key, dtype=bool)
        valid_studies = np.flatnonzero(valid_mask)
        n_valid = len(valid_studies)
        
        if n_valid - 2 <= num_pcs:
            for idx in variant_indices:
                results["Comments"][idx] = "SmallCohortCount"
            continue
        
        # Build design matrix X (same for all variants in this group)
        X = np.ones((n_valid, num_pcs + 1), dtype=np.float64)
        X[:, 1:] = pcs[valid_studies, :num_pcs]
        
        # Extract data for all variants in this batch
        n_batch = len(variant_indices)
        Y_batch = beta_data[variant_indices][:, valid_studies]  # (n_batch, n_valid)
        SE_batch = se_data[variant_indices][:, valid_studies]  # (n_batch, n_valid)
        W_batch = 1.0 / (SE_batch ** 2)  # (n_batch, n_valid)
        
        # Process batch using vectorized weighted least squares
        beta_hat_batch, se_beta_batch, TSS_batch, TSS0_batch, RSS_batch = _weighted_least_squares_batch(
            Y_batch, X, W_batch
        )
        
        if beta_hat_batch is None:
            for idx in variant_indices:
                results["Comments"][idx] = "collinear"
            continue
        
        # Vectorized processing of results for all variants in batch
        # Store coefficients (vectorized)
        result_beta_0[variant_indices] = beta_hat_batch[:, 0]
        result_se_0[variant_indices] = se_beta_batch[:, 0]
        
        # Store PC coefficients (vectorized)
        if num_pcs > 0:
            for j in range(num_pcs):
                result_beta_pc[j][variant_indices] = beta_hat_batch[:, j + 1]
                result_se_pc[j][variant_indices] = se_beta_batch[:, j + 1]
        
        # Vectorized statistics calculations
        TSS_RSS_batch = np.abs(TSS_batch - RSS_batch)  # (n_batch,)
        TSS0_RSS_batch = np.abs(TSS0_batch - RSS_batch)  # (n_batch,)
        chisq_assoc_batch = TSS0_RSS_batch  # (n_batch,)
        chisq_ancestry_batch = TSS_RSS_batch  # (n_batch,)
        chisq_residual_batch = np.abs(RSS_batch)  # (n_batch,)
        
        # Store association chi-squares for GCO (vectorized)
        valid_chisq_mask = chisq_assoc_batch > 0
        n_valid_chisq = valid_chisq_mask.sum()
        if n_valid_chisq > 0:
            end_idx = good_chisq_idx + n_valid_chisq
            good_chisq_assoc[good_chisq_idx:end_idx] = chisq_assoc_batch[valid_chisq_mask]
            good_chisq_idx = end_idx
        
        # Vectorized p-value calculations
        ndf_ancestry = num_pcs
        ndf_residual = n_valid - num_pcs_plus_1
        
        # Calculate p-values for ancestry heterogeneity (vectorized)
        p_ancestry_batch = np.full(n_batch, np.nan, dtype=np.float64)
        valid_ancestry_mask = chisq_ancestry_batch > 0
        if valid_ancestry_mask.any():
            p_ancestry_batch[valid_ancestry_mask] = chi2_sf(
                chisq_ancestry_batch[valid_ancestry_mask], ndf_ancestry
            )
        
        # Calculate p-values for residual heterogeneity (vectorized)
        p_residual_batch = np.full(n_batch, np.nan, dtype=np.float64)
        valid_residual_mask = (chisq_residual_batch > 0) & (ndf_residual > 0)
        if valid_residual_mask.any():
            p_residual_batch[valid_residual_mask] = chi2_sf(
                chisq_residual_batch[valid_residual_mask], ndf_residual
            )
        
        # Calculate lnBF (vectorized)
        lnBF_batch = (TSS0_RSS_batch - num_pcs_plus_1 * np.log(n_valid)) / 2.0
        
        # Store statistics (vectorized assignments)
        result_chisq_assoc[variant_indices] = chisq_assoc_batch
        result_ndf_assoc[variant_indices] = num_pcs_plus_1
        result_p_assoc[variant_indices] = np.nan  # Will be calculated later after GCO
        result_chisq_ancestry[variant_indices] = chisq_ancestry_batch
        result_ndf_ancestry[variant_indices] = ndf_ancestry
        result_p_ancestry[variant_indices] = p_ancestry_batch
        result_chisq_residual[variant_indices] = chisq_residual_batch
        result_ndf_residual[variant_indices] = ndf_residual
        result_p_residual[variant_indices] = p_residual_batch
        result_lnBF[variant_indices] = lnBF_batch
        
        # Calculate average EAF weighted by sample size (vectorized)
        if n_data is not None and eaf_data is not None:
            # Extract data for all variants in batch: (n_batch, n_valid)
            n_batch_data = n_data[variant_indices][:, valid_studies]  # (n_batch, n_valid)
            eaf_batch_data = eaf_data[variant_indices][:, valid_studies]  # (n_batch, n_valid)
            
            # Vectorized weighted EAF calculation
            total_n_batch = np.nansum(n_batch_data, axis=1)  # (n_batch,)
            weighted_eaf_numerator = np.nansum(n_batch_data * eaf_batch_data, axis=1)  # (n_batch,)
            
            # Only update where total_n > 0
            valid_n_mask = total_n_batch > 0
            if valid_n_mask.any():
                results["EAF"][variant_indices[valid_n_mask]] = (
                    weighted_eaf_numerator[valid_n_mask] / total_n_batch[valid_n_mask]
                )
                results["Nsample"][variant_indices[valid_n_mask]] = total_n_batch[valid_n_mask].astype(int)
        
        # Store Ncohort (vectorized)
        results["Ncohort"][variant_indices] = n_valid
        
        # Effects direction string (fully vectorized)
        # Extract beta data for all variants in batch: (n_batch, nstudy)
        beta_batch_all = beta_data[variant_indices, :]  # (n_batch, nstudy)
        
        # Create effects arrays for all variants at once: (n_batch, nstudy)
        effects_arrays = np.full((n_batch, nstudy), "?", dtype=object)
        
        # Vectorized effects direction calculation for valid studies
        beta_batch_valid = beta_batch_all[:, valid_studies]  # (n_batch, n_valid)
        effects_valid = np.select(
            [beta_batch_valid > 0, beta_batch_valid < 0, beta_batch_valid == 0],
            ["+", "-", "0"],
            default="?"
        )  # (n_batch, n_valid)
        
        # Assign valid study effects (vectorized)
        effects_arrays[:, valid_studies] = effects_valid
        
        # Join effects arrays into strings (vectorized using list comprehension for string operations)
        # Note: String operations are inherently sequential, but we've minimized the work
        # Convert numpy array slice to list for efficient string joining
        effects_strings = ["".join(effects_arrays[i, :].tolist()) for i in range(n_batch)]
        # Assign to list (results["Effects"] is a list, not numpy array)
        for i, idx in enumerate(variant_indices):
            results["Effects"][idx] = effects_strings[i]
    
    # Step 6: Calculate second GC correction (gco) if enabled
    # MR-MEGA has --gco option (separate from --gc) that applies a second genomic control correction
    # to the output association chi-squares
    # Filter good_chisq_assoc to remove NaN values (from pre-allocation)
    if good_chisq_idx > 0:
        good_chisq_assoc = good_chisq_assoc[:good_chisq_idx]
    else:
        good_chisq_assoc = np.array([], dtype=np.float64)
    
    gco_lambda = 1.0
    if use_gco and len(good_chisq_assoc) > 0:
        log.write("Calculating second genomic control (gco) lambda...", verbose=verbose)
        # Sort chi-squares and get median
        sorted_chisq = np.sort(good_chisq_assoc)
        n_chisq = len(sorted_chisq)
        if n_chisq % 2 == 0:
            median_chisq = (sorted_chisq[n_chisq // 2 - 1] + sorted_chisq[n_chisq // 2]) / 2.0
        else:
            median_chisq = sorted_chisq[(n_chisq - 1) // 2]
        
        # gco_lambda = median_chisq / median_of_chi2(df=num_pcs+1)
        # MR-MEGA uses invchisquaredistribution(num_pcs+1, 0.5) which is the median
        # chi2 is already imported at the top of the file
        median_expected = chi2.ppf(0.5, num_pcs + 1)  # median of chi2(df=num_pcs+1)
        gco_lambda = median_chisq / median_expected if median_expected > 0 else 1.0
        
        if gco_lambda > 1.0:
            log.write(f"GCO lambda: {gco_lambda:.4f}", verbose=verbose)
    
    if use_gco:
        log.write("Applying second GC correction (gco) and calculating final p-values...", verbose=verbose)
        chisq_assoc = results["chisq_association"]
        ndf_assoc = results["ndf_association"]
        valid_mask = ~(np.isnan(chisq_assoc) | (chisq_assoc <= 0) | np.isnan(ndf_assoc))
        
        chisq_assoc[valid_mask] /= gco_lambda
        
        # Vectorized p-value calculation (all valid variants have same ndf = num_pcs_plus_1)
        p_corrected = results["P-value_association"]
        if valid_mask.any():
            p_corrected[valid_mask] = chi2_sf(chisq_assoc[valid_mask], num_pcs_plus_1)
    else:
        log.write("Calculating p-values (gco correction disabled)...", verbose=verbose)
        chisq_assoc = results["chisq_association"]
        ndf_assoc = results["ndf_association"]
        valid_mask = ~(np.isnan(chisq_assoc) | (chisq_assoc <= 0) | np.isnan(ndf_assoc))
        
        # Vectorized p-value calculation (all valid variants have same ndf = num_pcs_plus_1)
        p_values = results["P-value_association"]
        if valid_mask.any():
            p_values[valid_mask] = chi2_sf(chisq_assoc[valid_mask], num_pcs_plus_1)
    
    # Create result DataFrame
    result_df = pd.DataFrame(results, index=data.index)
    
    # Add SNP info
    for col in ["SNPID", "CHR", "POS", "EA", "NEA"]:
        if col in data.columns:
            result_df[col] = data[col].values
    
    # Rename columns to match MR-MEGA output format
    column_mapping = {
        "SNPID": "MarkerName",
        "CHR": "Chromosome",
        "POS": "Position",
    }
    result_df = result_df.rename(columns=column_mapping)
    
    # Map to GWASLab format
    result_df["BETA"] = result_df["beta_0"]
    result_df["SE"] = result_df["se_0"]
    result_df["P"] = result_df["P-value_association"]
    result_df["Z"] = result_df["BETA"] / result_df["SE"]
    result_df["Q"] = result_df["chisq_residual_het"]
    result_df["DOF"] = result_df["ndf_residual_het"]
    result_df["P_HET"] = result_df["P-value_residual_het"]
    result_df["I2"] = np.where(
        result_df["Q"] > 0,
        np.maximum(0.0, (result_df["Q"] - result_df["DOF"]) / result_df["Q"]),
        0.0
    )
    
    # Create Sumstats object
    standard_cols = ["MarkerName", "Chromosome", "Position", "EA", "NEA", 
                     "BETA", "SE", "P", "Z", "EAF", "Nsample", "Q", "DOF", "P_HET", "I2"]
    other_cols = [col for col in result_df.columns if col not in standard_cols]
    
    # Rename back to GWASLab standard names
    result_df = result_df.rename(columns={"MarkerName": "SNPID", "Chromosome": "CHR", "Position": "POS"})
    
    result_sumstats = Sumstats(result_df, fmt="gwaslab", other=other_cols)
    
    log.write("MR-MEGA analysis completed successfully", verbose=verbose)
    
    return result_sumstats


def meta_regress_mrmega(
    sumstats_multi: 'SumstatsMulti',
    num_pcs: int = 2,
    use_genomic_control: bool = True,
    use_gco: bool = True,
    min_maf: float = 0.01,
    log: Optional[Log] = None,
    verbose: bool = True,
    selected_marker_snpids: Optional[list] = None
) -> Optional[Sumstats]:
    """
    Perform meta-regression using MR-MEGA algorithm (pure Python implementation).
    
    This function reimplements the MR-MEGA algorithm in pure Python without requiring
    the external executable. MR-MEGA uses an MDS-based approach:
    
    1. **MDS-based approach**: MR-MEGA calculates Multidimensional Scaling (MDS)
       axes from EAF (Effect Allele Frequency) correlation distance between
       cohorts. These MDS axes capture genetic diversity patterns.
    
    2. **Meta-regression model**: 
       Y = beta_0 + beta_1*PC1 + beta_2*PC2 + ... + error
       where:
       - Y = beta values from each cohort
       - PC1, PC2, ... = MDS principal components
       - Weights = 1/SE^2 (inverse variance, with GC correction if enabled)
    
    3. **Constraint**: Number of PCs must satisfy: npcs < ncohorts - 2
       (strict inequality from MR-MEGA source code)
    
    4. **Statistics**: MR-MEGA calculates:
       - Association test: chi-square for all coefficients (intercept + PCs)
       - Ancestry heterogeneity: chi-square for PC coefficients only
       - Residual heterogeneity: chi-square for residuals
    
    Parameters
    ----------
    sumstats_multi : SumstatsMulti
        SumstatsMulti object containing summary statistics from multiple studies.
    num_pcs : int, default 2
        Number of principal components to use. Must satisfy: num_pcs < nstudy - 2
    use_genomic_control : bool, default True
        Whether to apply genomic control correction to standard errors (--gc flag).
    use_gco : bool, default True
        Whether to apply second genomic control correction to output (--gco flag).
    min_maf : float, default 0.01
        Minimum minor allele frequency for marker selection (1%).
    log : Log, optional
        Log object for logging messages. If None, creates a new Log.
    verbose : bool, default True
        Whether to print progress messages.
    
    Returns
    -------
    Sumstats or None
        Sumstats object containing MR-MEGA results.
    """
    if log is None:
        log = Log()
    
    log.write(
        "Using MR-MEGA algorithm (pure Python implementation). "
        "MDS-based approach (PCs derived from EAF correlation).",
        verbose=verbose
    )
    
    # Call the Python implementation
    return meta_regress_mrmega_python(
        sumstats_multi=sumstats_multi,
        num_pcs=num_pcs,
        use_genomic_control=use_genomic_control,
        use_gco=use_gco,
        min_maf=min_maf,
        log=log,
        verbose=verbose
    )
