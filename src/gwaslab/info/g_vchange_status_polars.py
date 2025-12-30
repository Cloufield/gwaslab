from typing import Union, List
import pandas as pd
import polars as pl

def vchange_statusp(sumstats: pl.DataFrame, matched_index: pl.Expr, status: str, digit: int, before: List[str], after: List[str]) -> pl.DataFrame:
    dic={}
    for i in range(len(before)):
        dic[before[i]]=after[i]
    
    sumstats = sumstats.with_columns(pl.col(status).cast(pl.String).alias(status))

    if digit>1:
        sumstats = sumstats.with_columns(
            pl.when( matched_index )  
                .then(   pl.col(status).str.slice(0,digit-1) +   pl.col(status).str.slice(digit-1,1).str.replace_many(dic)  + pl.col(status).str.slice(digit))  
                .otherwise( pl.col(status) )
                .alias(status)
        )
    else:
        sumstats = sumstats.with_columns(
            pl.when( matched_index )  
                .then( pl.col(status).str.slice(0,1).str.replace_many(dic)  + pl.col(status).str.slice(digit)  )  
                .otherwise( pl.col(status) )
                .alias(status)
        )
    return sumstats

def copy_statusp(sumstats: pl.DataFrame, matched_index: pl.Expr, from_status: str, to_status: str, digit: int) -> pl.DataFrame:
    sumstats = sumstats.with_columns(pl.col(from_status).cast(pl.String).alias(from_status))
    sumstats = sumstats.with_columns(pl.col(to_status).cast(pl.String).alias(to_status))
    if digit>1:
        sumstats = sumstats.with_columns(
            pl.when( matched_index )  
                .then(   pl.col(from_status).str.slice(0,digit-1) +   pl.col(to_status).str.slice(digit-1,1)  + pl.col(from_status).str.slice(digit))  
                .otherwise( pl.col(to_status) )
                .alias(to_status)
        )
    else:
        sumstats = sumstats.with_columns(
            pl.when( matched_index )  
                .then( pl.col(from_status).str.slice(0,1)  + pl.col(to_status).str.slice(digit)  )  
                .otherwise( pl.col(to_status) )
                .alias(to_status)
        )
    return sumstats