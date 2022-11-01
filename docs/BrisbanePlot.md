#  Brisbane plot

gwaslab can create Brisbane plot (GWAS hits density plot).

See examples [here.](https://cloufield.github.io/gwaslab/visualization_brisbane/)

```
mysumstats = gl.Sumstats("height_lead.tsv",
             snpid="SNP",
             chrom="Chr",
             pos="BP_HG19",
             p="P-value")
             
mysumstats.plot_mqq(mode="b",anno="GENENAME",
                  build="19",anno_fixed_arm_length=2,
                  anno_args={"rotation":90},
                  marker_size=(30,30),sig_line_color="red")

```

- `mode="b"` : create Brisbane plot !
- `bwindowsizekb = 100` : windowsize in kb.

<img width=700 src="https://user-images.githubusercontent.com/40289485/197393168-e3e7076f-2801-4d66-9526-80778d44f3da.png">

Citation:Yengo, L., Vedantam, S., Marouli, E., Sidorenko, J., Bartell, E., Sakaue, S., ... & Lee, J. Y. (2022). A saturated map of common genetic variants associated with human height. Nature, 1-16.
