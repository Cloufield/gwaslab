# Extract lead variants

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

GWASLab will extract the lead variants from identified significant loci based on a sliding window (default window size: 500kb). (Details are described in )

Options:
- `windowsizekb` :specify the sliding window size in kb (default: 500)
- `sig_level` :specify the P value threshold (default: 5e-8).
- `xymt`: list of notation for chrX, chrY and chrMT.
- `anno`: `boolean` if annotate the lead variants with nearest gene names.
- `build` : `string` genome build version "19" or "38".

Return a dataframe of the lead variants

## Example

Sample sumstats: IS from pheweb.jp [https://pheweb.jp/pheno/IS](https://pheweb.jp/pheno/IS)

![image](https://user-images.githubusercontent.com/40289485/196447159-f9b41510-feb6-4ec8-adeb-a8f4c3683061.png)


```
mysumstats = gl.Sumstats("./hum0197.v3.BBJ.IS.v1/GWASsummary_IS_Japanese_SakaueKanai2020.auto.txt.gz",  fmt="saige")
mysumstats.get_lead()
```

![image](https://user-images.githubusercontent.com/40289485/196446293-c8f9c2d3-82c3-4122-bddd-184b43f2950f.png)

