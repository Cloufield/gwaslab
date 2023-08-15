# Regional plots

<img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/197126045-b1c55adf-3391-4c3d-b2f6-eaeac7c26024.png">

GWASLab provides functions for creating regional plots.

## 1. Usage
```
.plot_mqq(mode="r",
          region = None,
          ...
          ):
```

GWASLab regional plot function is based on plot_mqq().
Most options are largely the same as [Manhattan plot](https://cloufield.github.io/gwaslab/Visualization/).

|Option|DataType|Description|Default|
|-|-|-|-|
|`mode`|`r`|specify regional plot mode|-|
|`region`|`tuple`|a three elements tuple (chr, start, end); for example, (7,156538803,157538803)|-|
|`vcf_path`|`string`|path to LD reference in VCF format: if None, LD information will not be plotted.|`None`|
|`region_ref`|`boolean`|the SNPID or rsID (if SNPID is not available) for reference variant; if None, lead variants will be selected|`None`|
|`region_ref2`|`boolean`|the SNPID or rsID for the second reference variant|`None`|
|`region_grid`|`boolean`|If True, plot the grid line|`False`|
|`region_grid_line`|`dict`|parameters for the grid line|`{"linewidth": 2,"linestyle":"--"}`|
|`region_lead_grid`|`string`|If True, plot a line to show the reference variants|-|
|`region_lead_grid_line`|`string`|parameters for the line to show the reference variants|{"alpha":0.5,"linewidth" : 2,"linestyle":"--","color":"#FF0000"}|
|`region_ld_threshold`|`list`|LD r2 categories|`[0.2,0.4,0.6,0.8]`|
|`region_ld_colors`|`list`|LD r2 categories colors|`["#E4E4E4","#020080","#86CEF9","#24FF02","#FDA400","#FF0000","#FF0000"]`|
|`region_ld_colors1`|`list`|LD r2 categories colors for the first reference varaint (when region_ref2 is specified)|`["#E4E4E4","#F8CFCF","#F5A2A5","#F17474","#EB4445","#E51819","#E51819"]`|
|`region_ld_colors2`|`list`|LD r2 categories colors for the second reference varaint (when region_ref2 is specified)|`["#E4E4E4","#D8E2F2","#AFCBE3","#86B3D4","#5D98C4","#367EB7","#367EB7"]`|
|`region_hspace`|`float`|the space between the scatter plot and the gene track|`0.02`|
|`region_step`|`int`|number of X axis ticks|`21`|
|`region_recombination`|`boolean`||`True`|
|`tabix`|`string`|path to tabix; if None, GWASLab will search in environmental path; Note: if tabix is available, the speed is much faster!!!|`None`|
|`taf`|`list`|a five-element list; number of gene track lanes, offset for gene track, font_ratio, exon_ratio, text_offset|`[4,0,0.95,1,1]`|
|`build`|`19` or `38`|reference genome build; `99` for unknown|`99`|


!!! info "Calculation of LD r2"
    The calculation is based on [Rogers and Huff r implemented in scikit-alle](https://scikit-allel.readthedocs.io/en/stable/stats/ld.html). Variants in refernece vcf file should be biallelic format. Unphased data is acceptable. AF information is not needed. Variant ID is not required. Missing genotype is allowed.


## 2. Examples

### 2.1 Regional mqq plot
```
!wget -O t2d_bbj.txt.gz http://jenger.riken.jp/14/
```

```
mysumstats = gl.Sumstats("t2d_bbj.txt.gz",
             snpid="SNP",
             chrom="CHR",
             pos="POS",
             neaf="Frq",
             p="P")

mysumstats.plot_mqq(region=(7,156538803,157538803))
```
<img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/197125768-89c7fdd0-c80a-4db6-b8fc-9e970b39610e.png">


### 2.2 Regional plot without LD
```
mysumstats.plot_mqq(mode="r", region=(7,156538803,157538803),region_grid=True)
```
<img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/197126849-72439c54-19fe-4fb9-a795-e35f42a9236b.png">



### 2.3 Regional plot
```
mysumstats.plot_mqq(mode="r",region=(7,156538803,157538803),region_grid=True,
                    vcf_path="/home/yunye/mydata/d_disk/eas_1kg_af/EAS.chr7.split_norm_af.vcf.gz")
```
<img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/197126045-b1c55adf-3391-4c3d-b2f6-eaeac7c26024.png">

### 2.4 Two-reference-variant regional plot
```
mysumstats.plot_mqq(mode="r",region=(7,156538803,157538803),
                    region_grid=True,
                    region_ref ="7:156994964_C_T",
                    region_ref2="7:156793713_T_C",
                    vcf_path=gl.get_path("1kg_eas_hg19"))
```
<img width="600" alt="image" src="https://github.com/Cloufield/gwaslab/assets/40289485/930bed00-74ab-4870-947c-06bc11b932c7">
