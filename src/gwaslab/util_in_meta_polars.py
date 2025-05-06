
import pandas as pd
import numpy as np
from scipy.stats.distributions import chi2
from scipy.stats import norm
from gwaslab.g_Log import Log
from gwaslab.io_to_pickle import load_data_from_pickle
from gwaslab.g_Sumstats import Sumstats
import polars as pl
########################################################################################################################################################################################################################################################################################################################################################
########################################################################################################################################################################################################################################################################################################################################################
########################################################################################################################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
########################################################################################################################################################################################################################################################################################################################################################
########################################################################################################################################################################################################################################################################################################################################################

def meta_analyze_polars(sumstats_multi,
                       random_effects=False,
                       nstudy=1, 
                       log=Log()):
    log.write("Start to perform meta-analysis...")
    ###########################################################################
    log.write(" -Initiating result DataFrame...")
    
    sumstats_multi = sumstats_multi.with_row_index()

    sumstats_multi = sumstats_multi.with_columns(
        N=pl.lit(0),
        _BETAW_SUM = pl.lit(0.0),
        _BETA2W_SUM =pl.lit(0.0),
        _W_SUM =pl.lit(0.0),
        _W2_SUM=pl.lit(0.0),
        _EA_N=pl.lit(0.0),
        _NEA_N=pl.lit(0.0),
        DIRECTION =pl.lit(""),
        BETA=pl.lit(0.0),
        SE=pl.lit(0.0),
        DOF =pl.lit(-1),
        _R2=pl.lit(0.0),
    )

    for i in range(nstudy):
        n="N_{}".format(i+1)
        beta="BETA_{}".format(i+1)
        se="SE_{}".format(i+1)
        eaf="EAF_{}".format(i+1)
        sumstats_multi = sumstats_multi.with_columns(
        pl.when(pl.col(beta).is_null())  
        .then(pl.col("N"))  
        .otherwise(pl.col("N") + pl.col(n))
        .alias("N")
        ).with_columns(
        pl.when(pl.col(beta).is_null())  
        .then(pl.col("DOF") )  
        .otherwise(pl.col("DOF") + 1 )
        .alias("DOF")
        ).with_columns(   
        pl.when(pl.col(beta).is_null())  
        .then(pl.col("_BETA2W_SUM"))  
        .otherwise(pl.col("_BETA2W_SUM")   +   pl.col(beta)**2 *(1/(pl.col(se)**2))) 
        .alias("_BETA2W_SUM")
        ).with_columns(   
        pl.when(pl.col(beta).is_null())  
        .then( pl.col("_BETAW_SUM") )  
        .otherwise(pl.col("_BETAW_SUM") + pl.col(beta)*(1/(pl.col(se)**2)) )
        .alias("_BETAW_SUM")        
        ).with_columns(
        pl.when(pl.col(beta).is_null())  
        .then( pl.col("_W_SUM") )  
        .otherwise(pl.col("_W_SUM") + 1/(pl.col(se)**2))   
        .alias("_W_SUM")   
        ).with_columns(
        pl.when(pl.col(beta).is_null())  
        .then( pl.col("_W2_SUM") )  
        .otherwise(pl.col("_W2_SUM") +  pl.col("_W_SUM")**2)  
        .alias("_W2_SUM")   
        ).with_columns(
        pl.when(pl.col(beta).is_null())  
        .then( pl.col("_EA_N") )  
        .otherwise(pl.col("_EA_N") +  pl.col(n)*pl.col(eaf))  
        .alias("_EA_N")   
        ).with_columns(
        pl.when(pl.col(beta).is_null())  
        .then( pl.col("_NEA_N") )  
        .otherwise(pl.col("_EA_N") +  pl.col(n)*(1-pl.col(eaf)))  
        .alias("_NEA_N")   
        ).with_columns(
        pl.when(pl.col(beta) <0 )  
        .then( pl.col("DIRECTION") +"-" )  
        .otherwise( pl.col("DIRECTION") )  
        .alias("DIRECTION")   
        ).with_columns(
        pl.when(pl.col(beta) >0 )  
        .then( pl.col("DIRECTION") +"-" )  
        .otherwise( pl.col("DIRECTION") )  
        .alias("DIRECTION")   
        ).with_columns(
        pl.when(pl.col(beta) ==0 )  
        .then( pl.col("DIRECTION") +"0" )  
        .otherwise( pl.col("DIRECTION") )  
        .alias("DIRECTION")   
        ).with_columns(
        pl.when(pl.col(beta).is_null() )  
        .then( pl.col("DIRECTION") +"?" )  
        .otherwise( pl.col("DIRECTION") )  
        .alias("DIRECTION")   
        )
        
    sumstats_multi = sumstats_multi.with_columns(
        BETA = pl.col("_BETAW_SUM") / pl.col("_W_SUM"),
        EAF = pl.col("_EA_N") / (pl.col("_EA_N") + pl.col("_NEA_N")),
        SE = (1/pl.col("_W_SUM")).sqrt(),
        Q = pl.col("_BETA2W_SUM") - (pl.col("_BETAW_SUM")**2/pl.col("_W_SUM"))
    )

    sumstats_multi = sumstats_multi.with_columns(
        Z = pl.col("BETA") / pl.col("SE")
    ).with_columns(
        P = pl.col("Z").map_batches(lambda x: pl.Series(2*norm.sf(x.abs())))
    )

    if random_effects==True:
        sumstats_multi = sumstats_multi.with_columns(
            _R2 = (pl.col("Q") - pl.col("DOF")) / (pl.col("_W_SUM") - (pl.col("_W2_SUM")/pl.col("_W_SUM")))
        ).with_columns(
           pl.when(pl.col("_R2")<0 )  
            .then( pl.lit(0)  )  
            .otherwise( pl.col("_R2") )  
            .alias("_R2")   
        )
        
        sumstats_multi = sumstats_multi.with_columns(
            _BETAW_SUM_R = pl.lit(0.0),
            _W_SUM_R =pl.lit(0.0),
            BETA_RANDOM =pl.col("BETA"),
            SE_RANDOM =  pl.col("SE")
            )


        
        for i in range(nstudy):
            n="N_{}".format(i+1)
            beta="BETA_{}".format(i+1)
            se="SE_{}".format(i+1)
            eaf="EAF_{}".format(i+1)
            
            sumstats_multi = sumstats_multi.with_columns(
                pl.when( pl.col(beta).is_null() )  
                .then(pl.col("_BETAW_SUM_R"))  
                .otherwise(pl.col("_BETAW_SUM_R") + pl.col(beta)*(1/ (pl.col(se)**2 + pl.col("_R2")) ))
                .alias("_BETAW_SUM_R")
            ).with_columns(
                pl.when( pl.col(beta).is_null() )  
                .then(pl.col("_W_SUM_R") )  
                .otherwise( pl.col("_W_SUM_R")  + 1/(pl.col(se)**2 + pl.col("_R2") ) )
                .alias("_W_SUM_R")
            )

        sumstats_multi = sumstats_multi.with_columns(
            BETA_RANDOM = (pl.col("_BETAW_SUM_R") / pl.col("_W_SUM_R"))
        ).with_columns(
            SE_RANDOM = ( (1/pl.col("_W_SUM_R")).sqrt() )
        ).with_columns(
            Z_RANDOM = pl.col("BETA_RANDOM") /  pl.col("SE_RANDOM")
        ).with_columns(
            P_RANDOM = pl.col("Z_RANDOM").map_batches(lambda x: pl.Series(2*norm.sf(x.abs()))
        )
    )

    sumstats_multi = sumstats_multi.select(pl.all().exclude("^_.*$|\w+_[\d]+$")        )
    log.write("Finished performing meta-analysis.")
    return sumstats_multi