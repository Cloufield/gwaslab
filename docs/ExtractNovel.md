
# Checking if lead variants are novel or not

GWASLab can check if the lead variants of your summary statistics overlap with reported variants or not based on the physical distance.

## .get_novel()

```
sumstats.get_novel(
                    known,
                    efo,
                    only_novel=False,
                    windowsizekb_for_novel=1000,
                    windowsizekb=500,
                    sig_level=5e-8,
                    output_known=False)
```


GWASLab checks overlap with a local file of variants or records in GWASCatalog.

Required (either):

- `known` : `string`, path to the local file of reported variants

or 

- `efo` : `string`, EFO id for the target trait, which is used for querying the GWASCatalog.

## Options

| `.get_lead()` options    | DataType  | Description                                                                               | Default |
|--------------------------|-----------|-------------------------------------------------------------------------------------------|---------|
| `windowsizekb`           | `int`     | Specify the sliding window size in **kb**                                                 | `500`   |
| `sig_level`              | `float`   | Specify the P value threshold                                                             | `5e-8`  |
| `windowsizekb_for_novel` | `int`     | windowsize for determining if lead variants overlap with reported variants in GWASCatalog | `1000`  |
| `only_novel`             | `boolean` | output only novel variants                                                                | `False` |
| `output_known`           | `boolean` | additionally output the reported variants                                                 | `False` |
| `verbose`                | `boolean` | If True, print logs                                                                       | `True`  |


!!! info "EFO ID"
    You can find the efo id by simply searching in GWASCatalog. For example, the EFO ID for T2D can be obtained :
    
    <img width="700" alt="Screenshot 2023-02-03 at 13 36 16" src="https://user-images.githubusercontent.com/40289485/216513724-27f4e742-a03c-4346-a6bc-6b1002e22847.png">

## Example

!!! warning "Only works when your sumstats are based on GRCh38" 

!!! warning "Only associations with the EFO trait will be obtained. This does not include associations with child traits." 

!!! example 
    See [Lead and novel variants](https://cloufield.github.io/gwaslab/utility_get_lead_novel/)

## Reference

- Buniello, A., MacArthur, J. A. L., Cerezo, M., Harris, L. W., Hayhurst, J., Malangone, C., ... & Parkinson, H. (2019). The NHGRI-EBI GWAS Catalog of published genome-wide association studies, targeted arrays and summary statistics 2019. Nucleic acids research, 47(D1), D1005-D1012.
