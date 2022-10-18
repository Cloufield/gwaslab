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
```
mysumstats = gl.Sumstats(...)
mysumstats.basic_check()
mysumstats.get_lead()
```
