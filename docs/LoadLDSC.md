# Batch load LDSC log file

GWASLab provides a standalone function for batch loading ldsc log file into a `pd.DataFrame`. 

GWASLab uses regular expression to match the values and fill them into a dataframe.

## gl.read_ldsc()

```
gl.read_ldsc(filelist, mode="h2")
```

## Opions

- `filelist` : `list` .A list of paths to ldsc log files
- `mode` : `string`. `h2` or `rg`.  

## Example

!!! example ldsc - h2

    ```python
    #mode=h2
    paths = [
        '/home/he/work/GWASTutorial/08_LDSC/BBJ_HDLC.log',
        '/home/he/work/GWASTutorial/08_LDSC/BBJ_LDLC.log'
    ]
    
    df = gl.read_ldsc(ldsc_file_list, mode="h2")
    
    Loading file 1 :/home/he/work/GWASTutorial/08_LDSC/BBJ_HDLC.log ...
    Loading file 2 :/home/he/work/GWASTutorial/08_LDSC/BBJ_LDLC.log ...
    
    df
    Filename	h2_obs	h2_se	Lambda_gc	Mean_chi2	Intercept	Intercept_se	Ratio	Ratio_se
    BBJ_HDLC.log	0.1583	0.0281	1.1523	1.2843	1.0563	0.0114	0.1981	0.0402
    BBJ_LDLC.log	0.0743	0.0123	1.0833	1.1465	1.0296	0.0107	0.2019	0.0727
    ```
!!! example ldsc - rg    
    ```
    #mode=rg
    
    paths = [
        '/home/he/work/GWASTutorial/08_LDSC/BBJ_HDLC_LDLC.log',
    ]
    
    df = gl.read_ldsc(paths, mode="rg")
    Loading file 1 :/home/he/work/GWASTutorial/08_LDSC/BBJ_HDLC_LDLC.log ...
    
    df
    p1	p2	rg	se	z	p	h2_obs	h2_obs_se	h2_int	h2_int_se	gcov_int	gcov_int_se
    BBJ_HDLC.sumstats.gz	BBJ_LDLC.sumstats.gz	0.1601	0.1821	0.8794	0.3792	0.0543	0.0211	1.0583	0.0335	-0.0198	0.0121
    
    ```

## plot ldsc rg (under construction)

For genetic correlation, after loading, you can use `gl.plot_rg()` to plot a heat map to visualize the results. No extra manipulation needed.

```python
gl.plot_rg(myldscrg)
```
