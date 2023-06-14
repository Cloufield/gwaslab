Scatter & Distribution plot : allele frequency comparison


Available from v3.4.15

```
sumstats.check_af()

sumstats.plot_daf()
```

For check_af, see here.

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
    
    # check the difference in allele frequency with reference vcf
    mysumstats.check_af(ref_infer=gl.get_path("1kg_eas_hg19"), 
                        ref_alt_freq="AF",
                        n_cores=2)
    ```
    
    ```
    plot and get the outliers
    outliers = mysumstats.plot_daf(threshold=0.12, 
                                    save="af_correlation.png",
                                    save_args={"dpi":300})
    ```
    
    ```
    outliers[1]
    ```
