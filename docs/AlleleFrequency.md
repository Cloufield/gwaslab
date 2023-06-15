Scatter & Distribution plot : allele frequency comparison

Available from v3.4.15

```
#check the difference between the EAF in the sumstats and the allele frequency in VCF files
sumstats.check_af()

#allele frequnecy correlation plot
sumstats.plot_daf()
```

You need to run 'check_af()' first before plotting. For check_af(), see [here](https://cloufield.github.io/gwaslab/Harmonization/#check-the-difference-in-allele-frequency).

Options for `plot_daf`:
`threshold`: `float`, the threshold used to determine outliers.


!!! example
    ```
    mysumstats = gl.Sumstats("t2d_bbj.txt.gz",
                 snpid="SNP",
                 chrom="CHR",
                 pos="POS",
                 ea="ALT",
                 nea="REF",
                 neaf="Frq",
                 beta="BETA",
                 se="SE",
                 p="P",
                 direction="Dir",
                 n="N",nrows=10000)
    
    # harmonize
    mysumstats.harmonize(basic_check = True, 
                         ref_seq=gl.get_path("ucsc_genome_hg19"))
    ```
    <img width=600 src="https://github.com/Cloufield/gwaslab/assets/40289485/dc3716d6-1b85-42dd-99ca-d0a38a3ad6a1">

    ```
    # check the difference in allele frequency with reference vcf
    mysumstats.check_af(ref_infer=gl.get_path("1kg_eas_hg19"), 
                        ref_alt_freq="AF",
                        n_cores=2)
    ```
    <img width=600 src="https://github.com/Cloufield/gwaslab/assets/40289485/8cdc53b4-f661-40d6-a61b-ab5f23af6a4c">

    ```
    plot and get the outliers
    outliers = mysumstats.plot_daf(threshold=0.12, 
                                    save="af_correlation.png",
                                    save_args={"dpi":300})
    ```
    <img width=600 src="https://github.com/Cloufield/gwaslab/assets/40289485/0e9c92d4-3bea-4734-8ba6-76f7fc6e329b">

    ```
    outliers[1]
    ```
    <img width=600 src="https://github.com/Cloufield/gwaslab/assets/40289485/c99d91cb-f1f8-4412-bfbf-960457fc9d0e">

