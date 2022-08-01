# Utilities

- Data conversion

- Extract lead variants

- Format sumstats for commonly used tools

- Load LDSC results into pd.DataFrame

- Observed-scale heritability to liability-scale heritablity conversion

## Extract lead variants

```python
sumstats.get_lead(windowsizekb=500, sig_level=5e-8)
```

GWASLab will extract the lead variants from identified significant loci based on a sliding window (default window size: 500kb). (Details are described in )

`windowsizekb` :specify the sliding window size in kb (default: 500)

`sig_level` :specify the P value threshold (default: 5e-8)

Return a dataframe of the lead variants

## 

## Formatting sumstats

```
sumstats.format()
```

Format the sumstats to the formats that were accepted by commonly used tools including LDSC,  MAGMA, FUMA , METAL, (bed ,vcf,...) 

`path`: output file path

`format`: Currently, GWASlab support ldsc,   (fuma, metal,bed, vcf ...)  

`extract`: a list of SNPIDs. 

`exclude`: a list of SNPIDs.

`exclude_hla`: True or False. If True, exclude HLA region when exporting sumstats.

`hapmap3` : True or False. If True,, only exporting Hapmap3 SNPs.

`to_csvargs`: arguments for `pandas.to_csv() `function

## 

## Batch load LDSC log file

Simply batch load LDSc data into a pandas DataFrame for other manipulation.

GWASLab uses regular expression to match the values and fill them into a dataframe.

```gl.read_ldsc()```
gl.read_ldsc(filelist, mode="h2")

```
`filelist` : a list of paths to ldsc log files

`mode` : h2 or rg 

```python
#mode=h2
ldsc_file_list=["file1."]
myldscrg = gl.read_ldsc(ldsc_file_list, mode=h2)


#mode=rg
ldsc_file_list=["file1."]
myldscrg = gl.read_ldsc(mode=h2)
```

For genetic correlation, after loading, you can use `gl.plot_rg()`  to plot a heat map to visualize the results. No extra manipulation needed.

```python
gl.plot_rg(myldscrg)
```

## Heritabilty conversion (Observed-scale -> Liability-scale)

```
gl.h2_obs_to_liab(h2_obs, P, K, se_obs=None)
```

`h2_obs` : float. Heritability on the observed scale in an ascertained sample. 
`P `: float in (0,1).  Prevalence of the phenotype in the sample. 
`K `: float in (0,1) . Prevalence of the phenotype in the population. 
`se_obs` : float. se of h2_obs.

Adopted from LDSC.  

Reference: Estimating Missing Heritability for Disease from Genome-wide Association Studies  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3059431/ 
