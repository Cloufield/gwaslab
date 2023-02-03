# Extract lead variants

GWASLab will extract the lead variants from identified significant loci based on a sliding window.

!!! quote
    GWASLab basically adopted the definition for novel loci from Global Biobank Meta-analysis Initiative flagship paper. 
    
    **"We defined genome-wide significant loci by iteratively spanning the ±500 kb region around the most significant variant and merging overlapping regions until no genome-wide significant variants were detected within ±1 Mb."** 
    
    (Details are described in Zhou, W., Kanai, M., Wu, K. H. H., Rasheed, H., Tsuo, K., Hirbo, J. B., ... & Study, C. O. H. (2022). Global Biobank Meta-analysis Initiative: Powering genetic discovery across human disease. Cell Genomics, 2(10), 100192. )

    GWASlab currently iteratively extends ± `windowsizekb` kb region around the most significant variant and merges overlapping regions until no genome-wide significant variants were detected within ± `windowsizekb`. (slightly different from the GBMI paper. When `windowsizekb=1000`, it is equvalent to GBMI's definition.)

# Usage

```python
mysumstats.get_lead(
           windowsizekb=500,
           sig_level=5e-8,
           xymt=["X","Y","MT"],
           anno=False,
           build="19",
           source="ensembl",
           verbose=True)
```

## Options

- `windowsizekb` : `int`. Specify the sliding window size in **kb** (default: 500)
- `sig_level` : `float`.Specify the P value threshold (default: 5e-8).
- `xymt`: `list`. A list of notation for chrX, chrY and chrMT.
- `anno`: `boolean`. if annotate the lead variants with nearest gene names.
- `build` : `string`. genome build version "19" or "38".

Return a dataframe of the lead variants

## Example

!!! example
    Sample sumstats: IS from pheweb.jp [https://pheweb.jp/pheno/IS](https://pheweb.jp/pheno/IS)
    ![image](https://user-images.githubusercontent.com/40289485/196447159-f9b41510-feb6-4ec8-adeb-a8f4c3683061.png)
    ```
    mysumstats = gl.Sumstats("./hum0197.v3.BBJ.IS.v1/GWASsummary_IS_Japanese_SakaueKanai2020.auto.txt.gz",  fmt="saige")
    mysumstats.get_lead()
    ```
    ![image](https://user-images.githubusercontent.com/40289485/196446293-c8f9c2d3-82c3-4122-bddd-184b43f2950f.png)

