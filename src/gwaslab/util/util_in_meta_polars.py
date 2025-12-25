
import pandas as pd
import numpy as np
from scipy.stats.distributions import chi2
from scipy.stats import norm
from gwaslab.info.g_Log import Log
from gwaslab.io.io_to_pickle import load_data_from_pickle
from gwaslab.g_Sumstats import Sumstats
import polars as pl

############################################################################################################################################################################

def meta_analyze_polars(sumstats_multi,
                       random_effects=False,
                       nstudy=1, 
                       log=Log()):
    """
    Perform meta-analysis using polars DataFrame.
    
    Parameters
    ----------
    sumstats_multi : pl.DataFrame
        Polars DataFrame with columns BETA_1, SE_1, N_1, EAF_1, etc. for each study
    random_effects : bool, optional
        Whether to compute random effects model, by default False
    nstudy : int, optional
        Number of studies, by default 1
    log : Log, optional
        Log object for logging, by default Log()
    
    Returns
    -------
    pl.DataFrame
        Polars DataFrame with meta-analysis results
    """
    log.write("Start to perform meta-analysis...")
    log.write(" -Initiating result DataFrame...")
    
    # Add row index for tracking
    sumstats_multi = sumstats_multi.with_row_index()

    # Initialize accumulation columns
    sumstats_multi = sumstats_multi.with_columns([
        pl.lit(0).alias("N"),
        pl.lit(0.0).alias("_BETAW_SUM"),
        pl.lit(0.0).alias("_BETA2W_SUM"),
        pl.lit(0.0).alias("_W_SUM"),
        pl.lit(0.0).alias("_W2_SUM"),
        pl.lit(0.0).alias("_EA_N"),
        pl.lit(0.0).alias("_NEA_N"),
        pl.lit("").alias("DIRECTION"),
        pl.lit(0.0).alias("BETA"),
        pl.lit(0.0).alias("SE"),
        pl.lit(-1).alias("DOF"),
        pl.lit(0.0).alias("_R2"),
    ])

    # Iterate through each study to accumulate statistics
    for i in range(nstudy):
        n = f"N_{i+1}"
        beta = f"BETA_{i+1}"
        se = f"SE_{i+1}"
        eaf = f"EAF_{i+1}"
        
        # Calculate weight (inverse variance)
        weight = 1 / (pl.col(se) ** 2)
        beta_is_null = pl.col(beta).is_null()
        
        # Accumulate statistics in a single with_columns call for better performance
        sumstats_multi = sumstats_multi.with_columns([
            # N: sum of sample sizes
            pl.when(beta_is_null)
            .then(pl.col("N"))
            .otherwise(pl.col("N") + pl.col(n))
            .alias("N"),
            
            # DOF: degrees of freedom (number of studies with data)
            pl.when(beta_is_null)
            .then(pl.col("DOF"))
            .otherwise(pl.col("DOF") + 1)
            .alias("DOF"),
            
            # _BETA2W_SUM: sum of beta^2 * weight
            pl.when(beta_is_null)
            .then(pl.col("_BETA2W_SUM"))
            .otherwise(pl.col("_BETA2W_SUM") + pl.col(beta) ** 2 * weight)
            .alias("_BETA2W_SUM"),
            
            # _BETAW_SUM: sum of beta * weight
            pl.when(beta_is_null)
            .then(pl.col("_BETAW_SUM"))
            .otherwise(pl.col("_BETAW_SUM") + pl.col(beta) * weight)
            .alias("_BETAW_SUM"),
            
            # _W_SUM: sum of weights
            pl.when(beta_is_null)
            .then(pl.col("_W_SUM"))
            .otherwise(pl.col("_W_SUM") + weight)
            .alias("_W_SUM"),
            
            # _W2_SUM: sum of squared cumulative weights (for random effects)
            # Note: This accumulates (sum(w))^2 at each step, used in tau^2 calculation
            pl.when(beta_is_null)
            .then(pl.col("_W2_SUM"))
            .otherwise(pl.col("_W2_SUM") + pl.col("_W_SUM") ** 2)
            .alias("_W2_SUM"),
            
            # _EA_N: sum of N * EAF (effect allele count)
            pl.when(beta_is_null)
            .then(pl.col("_EA_N"))
            .otherwise(pl.col("_EA_N") + pl.col(n) * pl.col(eaf))
            .alias("_EA_N"),
            
            # _NEA_N: sum of N * (1 - EAF) (non-effect allele count)
            pl.when(beta_is_null)
            .then(pl.col("_NEA_N"))
            .otherwise(pl.col("_NEA_N") + pl.col(n) * (1 - pl.col(eaf)))
            .alias("_NEA_N"),
            
            # DIRECTION: build direction string based on beta sign
            # Order matters: check null first, then < 0, > 0, == 0
            pl.when(beta_is_null)
            .then(pl.col("DIRECTION") + "?")
            .when(pl.col(beta) < 0)
            .then(pl.col("DIRECTION") + "-")
            .when(pl.col(beta) > 0)
            .then(pl.col("DIRECTION") + "+")
            .when(pl.col(beta) == 0)
            .then(pl.col("DIRECTION") + "0")
            .otherwise(pl.col("DIRECTION"))
            .alias("DIRECTION"),
        ])
        
    # Calculate fixed-effects meta-analysis statistics
    sumstats_multi = sumstats_multi.with_columns([
        # BETA: weighted average
        (pl.col("_BETAW_SUM") / pl.col("_W_SUM")).alias("BETA"),
        
        # EAF: weighted average allele frequency
        (pl.col("_EA_N") / (pl.col("_EA_N") + pl.col("_NEA_N"))).alias("EAF"),
        
        # SE: standard error of weighted average
        (1 / pl.col("_W_SUM")).sqrt().alias("SE"),
        
        # Q: Cochran's Q statistic for heterogeneity
        (pl.col("_BETA2W_SUM") - (pl.col("_BETAW_SUM") ** 2 / pl.col("_W_SUM"))).alias("Q"),
    ])

    # Calculate Z-score and P-value
    sumstats_multi = sumstats_multi.with_columns([
        (pl.col("BETA") / pl.col("SE")).alias("Z"),
    ]).with_columns([
        pl.col("Z").map_batches(lambda x: pl.Series(2 * norm.sf(x.abs()))).alias("P"),
    ])
    
    # Calculate heterogeneity statistics (P_HET and I2)
    # I2: I-squared statistic for heterogeneity
    sumstats_multi = sumstats_multi.with_columns([
        ((pl.col("Q") - pl.col("DOF")) / pl.col("Q")).alias("I2"),
    ]).with_columns([
        pl.when(pl.col("I2") < 0)
        .then(pl.lit(0.0))
        .otherwise(pl.col("I2"))
        .alias("I2"),
    ])
    
    # P_HET: P-value for heterogeneity test (chi-square test)
    # Group by DOF and calculate P_HET for each group
    # Initialize P_HET column first
    sumstats_multi = sumstats_multi.with_columns([
        pl.lit(None).cast(pl.Float64).alias("P_HET"),
    ])
    
    dof_values = sumstats_multi["DOF"].unique().to_list()
    for dof in dof_values:
        if dof is not None and dof > 0:
            dof_mask = pl.col("DOF") == dof
            sumstats_multi = sumstats_multi.with_columns([
                pl.when(dof_mask)
                .then(pl.col("Q").map_batches(lambda x: pl.Series(chi2.sf(x, dof))))
                .otherwise(pl.col("P_HET"))
                .alias("P_HET"),
            ])

    # Random effects model
    if random_effects:
        log.write(" -Computing random effects model...")
        
        # Calculate tau^2 (between-study variance)
        sumstats_multi = sumstats_multi.with_columns([
            ((pl.col("Q") - pl.col("DOF")) / (pl.col("_W_SUM") - (pl.col("_W2_SUM") / pl.col("_W_SUM")))).alias("_R2"),
        ]).with_columns([
            pl.when(pl.col("_R2") < 0)
            .then(pl.lit(0.0))
            .otherwise(pl.col("_R2"))
            .alias("_R2"),
        ])
        
        # Initialize random effects accumulators
        sumstats_multi = sumstats_multi.with_columns([
            pl.lit(0.0).alias("_BETAW_SUM_R"),
            pl.lit(0.0).alias("_W_SUM_R"),
            pl.col("BETA").alias("BETA_RANDOM"),
            pl.col("SE").alias("SE_RANDOM"),
        ])

        # Recalculate weights with tau^2 added to variance
        for i in range(nstudy):
            beta = f"BETA_{i+1}"
            se = f"SE_{i+1}"
            beta_is_null = pl.col(beta).is_null()
            
            # Weight for random effects: 1 / (SE^2 + tau^2)
            weight_r = 1 / (pl.col(se) ** 2 + pl.col("_R2"))
            
            sumstats_multi = sumstats_multi.with_columns([
                pl.when(beta_is_null)
                .then(pl.col("_BETAW_SUM_R"))
                .otherwise(pl.col("_BETAW_SUM_R") + pl.col(beta) * weight_r)
                .alias("_BETAW_SUM_R"),
                
                pl.when(beta_is_null)
                .then(pl.col("_W_SUM_R"))
                .otherwise(pl.col("_W_SUM_R") + weight_r)
                .alias("_W_SUM_R"),
            ])

        # Calculate random effects statistics
        sumstats_multi = sumstats_multi.with_columns([
            (pl.col("_BETAW_SUM_R") / pl.col("_W_SUM_R")).alias("BETA_RANDOM"),
            ((1 / pl.col("_W_SUM_R")).sqrt()).alias("SE_RANDOM"),
        ]).with_columns([
            (pl.col("BETA_RANDOM") / pl.col("SE_RANDOM")).alias("Z_RANDOM"),
        ]).with_columns([
            pl.col("Z_RANDOM").map_batches(lambda x: pl.Series(2 * norm.sf(x.abs()))).alias("P_RANDOM"),
        ])

    # Remove temporary columns and study-specific columns
    # Keep only meta-analysis results and info columns
    sumstats_multi = sumstats_multi.select(
        pl.all().exclude(r"^_.*$|\w+_[\d]+$")
    )
    
    log.write("Finished performing meta-analysis.")
    return sumstats_multi
