# Liftover

Gwaslab can directly liftover the positions of variants in sumstats.

See examples [here.](https://cloufield.github.io/gwaslab/harmonization_liftover/)

```
mysumstats.liftover(n_cores=3, 
                    from_build="19", 
                    to_build="38",
                    remove=True)
```

## Options

- `n_cores` : number of threads to use for liftover
- `from_build` :  from which genome build
- `to_build` : to which genome build
- `remove` : remove unmapped variants

Perform liftover for POS (based on [liftover](https://github.com/jeremymcrae/liftover) GitHub - jeremymcrae/liftover: liftover for python, made fast with cython)

## example
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
             p="P",nrows=500000)
             
mysumstats.basic_check()

mysumstats.liftover(n_cores=3, from_build="19", to_build="38")

```

