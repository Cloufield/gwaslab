# Stacked Manhattan and regional plot

## Load gwaslab and sumstats

```python
import gwaslab as gl
```

```python
gl.show_version()
```

**stdout:**
```
2025/12/25 23:54:24 GWASLab v4.0.0 https://cloufield.github.io/gwaslab/
2025/12/25 23:54:24 (C) 2022-2025, Yunye He, Kamatani Lab, GPL-3.0 license, gwaslab@gmail.com
2025/12/25 23:54:24 Python version: 3.12.0 | packaged by conda-forge | (main, Oct  3 2023, 08:43:22) [GCC 12.3.0]
```

```python
gl1 = gl.Sumstats("../0_sample_data/bbj_bmi_female.txt.gz",beta="BETA",se="SE",snpid="SNP",chrom="CHR",ea="ALT",nea="REF",pos="POS",p="P", build="19",sep="\t", verbose=False)
gl2 = gl.Sumstats("../0_sample_data/bbj_bmi_male.txt.gz",snpid="SNP",chrom="CHR",ea="ALT",nea="REF",pos="POS",p="P", build="19",sep="\t", verbose=False)

gl1.basic_check(verbose=False)
gl2.basic_check(verbose=False)
```

| SNPID | CHR | POS | EA | NEA | STATUS | P |
| --- | --- | --- | --- | --- | --- | --- |
| rs143225517 | 1 | 751756 | C | T | 1980099 | 0.5883 |
| rs3094315 | 1 | 752566 | A | G | 1980099 | 0.5889 |
| rs61770173 | 1 | 753405 | A | C | 1980099 | 0.5878 |
| rs117086422 | 1 | 845635 | T | C | 1980099 | 0.5565 |
| rs28612348 | 1 | 846078 | T | C | 1980099 | 0.5659 |
| ... | ... | ... | ... | ... | ... | ... |
| rs41281537 | 22 | 51171667 | A | G | 1980099 | 0.8453 |
| rs756638 | 22 | 51171693 | A | G | 1980099 | 0.7726 |
| rs5770824 | 22 | 51172460 | C | T | 1980099 | 0.8374 |
| rs3810648 | 22 | 51175626 | G | A | 1980099 | 0.2884 |
| rs2285395 | 22 | 51178090 | A | G | 1980099 | 0.3223 |

*[11 rows x 7 columns]*

```python
gl.plot_stacked_mqq(objects=[gl1,gl2],
                    vcfs=[gl.get_path("1kg_eas_hg19"), gl.get_path("1kg_eas_hg19")],
                    region_legend_marker=False,
                    region=(19,46214297 - 100000, 46214297 + 100000),
                    region_ref =["rs35560038","rs1055220","rs4802274"],
                    build="19",
                    anno="SNPID",
                    fontsize=15,
                    mode="r",
                    region_chromatin_files=["../0_sample_data/E098_15_coreMarks_mnemonics.bed.gz"],
                    region_chromatin_labels=["E098_15"],
                    cbar_fontsize=11,
                    marker_size=(100,100),
                    titles=["Male","Female"],
                    title_kwargs={"size":20},
                    anno_kwargs={"rotation":0,"fontsize":10}, 
                    verbose=False, check=False)
```

**stdout:**
```
2025/12/25 23:55:37  -POS data type is already numeric. Skipping conversion...
2025/12/25 23:55:37  -CHR data type is already numeric. Skipping conversion...
2025/12/25 23:55:47  -POS data type is already numeric. Skipping conversion...
2025/12/25 23:55:47  -CHR data type is already numeric. Skipping conversion...
```

```
(<Figure size 3200x3200 with 7 Axes>,
 <gwaslab.info.g_Log.Log at 0x7fd21de3fb30>)
```

![Output image](images/notebooks/visualization_stacked_mqq_img_0.png)

## Stacked Manhattan plot

```python
gl.plot_stacked_mqq(objects=[gl1,gl2],
                    mode="m",
                    skip=3,
                    titles=["Female - BMI","Male - BMI"],
                    check=False,
                    verbose=False)
```

**stdout:**
```
2025/12/25 23:55:52 #WARNING! Genomic coordinates version is unknown.
2025/12/25 23:55:54  -POS data type is already numeric. Skipping conversion...
2025/12/25 23:55:54  -CHR data type is already numeric. Skipping conversion...
2025/12/25 23:55:57 #WARNING! Genomic coordinates version is unknown.
2025/12/25 23:55:59  -POS data type is already numeric. Skipping conversion...
2025/12/25 23:55:59  -CHR data type is already numeric. Skipping conversion...
```

```
(<Figure size 2000x1600 with 2 Axes>,
 <gwaslab.info.g_Log.Log at 0x7fd21de3fb30>)
```

![Output image](images/notebooks/visualization_stacked_mqq_img_1.png)

## Multiple stacked plots

```python
#gl.download_ref("1kg_eas_hg19")
```

```python
gl.plot_stacked_mqq(objects=[gl1,gl1,gl1,gl1],
                    mode="m",
                    skip=3,
                    titles=["BMI 1","BMI 2", "BMI 3", "BMI 4"],
                    check=False,
                    verbose=False)
```

**stdout:**
```
2025/12/25 23:56:03 #WARNING! Genomic coordinates version is unknown.
2025/12/25 23:56:05  -POS data type is already numeric. Skipping conversion...
2025/12/25 23:56:05  -CHR data type is already numeric. Skipping conversion...
2025/12/25 23:56:07 #WARNING! Genomic coordinates version is unknown.
2025/12/25 23:56:09  -POS data type is already numeric. Skipping conversion...
2025/12/25 23:56:09  -CHR data type is already numeric. Skipping conversion...
2025/12/25 23:56:12 #WARNING! Genomic coordinates version is unknown.
2025/12/25 23:56:14  -POS data type is already numeric. Skipping conversion...
2025/12/25 23:56:14  -CHR data type is already numeric. Skipping conversion...
2025/12/25 23:56:16 #WARNING! Genomic coordinates version is unknown.
2025/12/25 23:56:18  -POS data type is already numeric. Skipping conversion...
2025/12/25 23:56:18  -CHR data type is already numeric. Skipping conversion...
```

```
(<Figure size 2000x3200 with 4 Axes>,
 <gwaslab.info.g_Log.Log at 0x7fd21de3fb30>)
```

![Output image](images/notebooks/visualization_stacked_mqq_img_2.png)

## Add title and annotation

```python
gl.plot_stacked_mqq(objects=[gl1,gl2],
                    mode="m",
                    skip=3,
                    anno=True,
                    anno_style="right",
                    titles=["Female - BMI","Male - BMI"],
                    region_hspace=0.4,   # increase the space between each subplot
                    xpadl=0.05,          # increase left padding for X axis
                    title_pos=[0, 1.7],  # adjust title position
                    check=False,
                    verbose=False)
```

**stdout:**
```
2025/12/25 23:56:22 #WARNING! Genomic coordinates version is unknown.
2025/12/25 23:56:24  -POS data type is already numeric. Skipping conversion...
2025/12/25 23:56:24  -CHR data type is already numeric. Skipping conversion...
2025/12/25 23:56:27 #WARNING! Genomic coordinates version is unknown.
2025/12/25 23:56:28  -POS data type is already numeric. Skipping conversion...
2025/12/25 23:56:28  -CHR data type is already numeric. Skipping conversion...
```

```
(<Figure size 2000x1600 with 2 Axes>,
 <gwaslab.info.g_Log.Log at 0x7fd21de3fb30>)
```

![Output image](images/notebooks/visualization_stacked_mqq_img_3.png)

## Customization

You can customize each plot using options for plot_mqq() by adding a number to the end of the option like `highlight1` means `highlight` option for the first panel. if no number suffix was added, the option will be applied to all panels. 


```python
gl.plot_stacked_mqq(objects=[gl1,gl2],
                    mode="m",
                    skip=3,
                    anno=True,
                    anno_style="expand",
                    titles=["Female - BMI","Male - BMI"],
                    highlight1=["rs183975233"],
                    highlight2=["rs35560038"],
                    anno_set1=["rs183975233"],
                    anno_set2=["rs35560038"],
                    region_hspace=0.4,   # increase the space between each subplot
                    xpadl=0.05,          # increase left padding for X axis
                    title_pos=[0, 1.3],  # adjust title position
                    check=False,
                    verbose=False)
```

**stdout:**
```
2025/12/25 23:56:33 #WARNING! Genomic coordinates version is unknown.
2025/12/25 23:56:35  -POS data type is already numeric. Skipping conversion...
2025/12/25 23:56:35  -CHR data type is already numeric. Skipping conversion...
2025/12/25 23:56:40 #WARNING! Genomic coordinates version is unknown.
2025/12/25 23:56:41  -POS data type is already numeric. Skipping conversion...
2025/12/25 23:56:41  -CHR data type is already numeric. Skipping conversion...
```

```
(<Figure size 2000x1600 with 2 Axes>,
 <gwaslab.info.g_Log.Log at 0x7fd21de3fb30>)
```

![Output image](images/notebooks/visualization_stacked_mqq_img_4.png)

## Stacked regional plot

```python
gl2.filter_flanking_by_chrpos((19,46214297)).filter_value("P <1e-20").data
```

**stdout:**
```
2025/12/25 23:56:48 Start to extract variants in the flanking regions using CHR and POS ...(v4.0.0)
2025/12/25 23:56:48  -Current Dataframe shape : 5961600 x 7 ; Memory usage: 255.84 MB
2025/12/25 23:56:48  - Central positions: (19, 46214297)
2025/12/25 23:56:48  - Flanking windowsize in kb: 500
2025/12/25 23:56:48  - Variants in flanking region 19:45714297-46714297 : 2137
2025/12/25 23:56:48  - Extracted 2137 variants in the regions.
2025/12/25 23:56:48  -Filtered out variants: 5959463
2025/12/25 23:56:48  -Current Dataframe shape : 2137 x 7 ; Memory usage: 0.11 MB
2025/12/25 23:56:48 Finished extracting variants in the flanking regions.
2025/12/25 23:56:48 Start to filter variants by condition... ...(v4.0.0)
2025/12/25 23:56:48  -Current Dataframe shape : 2137 x 7 ; Memory usage: 0.11 MB
2025/12/25 23:56:48  -Expression: P <1e-20
2025/12/25 23:56:49  -Filtered out variants: 2131
2025/12/25 23:56:49  -Current Dataframe shape : 6 x 7 ; Memory usage: 0.00 MB
2025/12/25 23:56:49 Finished filtering variants.
```

| SNPID | CHR | POS | EA | NEA | STATUS | P |
| --- | --- | --- | --- | --- | --- | --- |
| rs35560038 | 19 | 46175046 | T | A | 1980099 | 1.643000e-24 |
| rs34089191 | 19 | 46176723 | C | G | 1980099 | 5.961000e-21 |
| rs55669001 | 19 | 46177235 | C | T | 1980099 | 9.950000e-21 |
| rs9676912 | 19 | 46212182 | T | C | 1980099 | 5.247000e-22 |
| rs8101428 | 19 | 46212952 | C | T | 1980099 | 1.094000e-22 |
| rs1055220 | 19 | 46214297 | T | C | 1980099 | 5.731000e-24 |

*[6 rows x 7 columns]*

## stacked regional plots with single reference SNP

```python
gl.plot_stacked_mqq(objects=[gl1,gl2],
                    vcfs=[gl.get_path("1kg_eas_hg19"), gl.get_path("1kg_eas_hg19")],
                    region=(19,46214297 - 300000, 46214297 + 300000),
                    build="19",
                    mode="r",
                    anno="SNPID",
                    anno_style="tight",
                    titles=["Male","Female"],
                    title_kwargs={"size":20},
                    anno_kwargs={"rotation":0}, verbose=False, check=False)
```

**stdout:**
```
2025/12/25 23:56:52  -POS data type is already numeric. Skipping conversion...
2025/12/25 23:56:52  -CHR data type is already numeric. Skipping conversion...
2025/12/25 23:57:05  -POS data type is already numeric. Skipping conversion...
2025/12/25 23:57:05  -CHR data type is already numeric. Skipping conversion...
```

```
(<Figure size 3200x2400 with 6 Axes>,
 <gwaslab.info.g_Log.Log at 0x7fd21de3fb30>)
```

![Output image](images/notebooks/visualization_stacked_mqq_img_5.png)

## stacked regional plots with multiple reference SNPs

```python
gl.plot_stacked_mqq(objects=[gl1,gl2],
                    vcfs=[gl.get_path("1kg_eas_hg19"), gl.get_path("1kg_eas_hg19")],
                    region=(19,46214297 - 100000, 46214297 + 100000),
                    region_ref =["rs35560038","rs1055220","rs4802274"],
                    build="19",
                    anno="SNPID",
                    anno_set =["rs35560038","rs1055220","rs4802274"],
                    mode="r",
                    anno_style="tight",anno_adjust=True,
                    cbar_fontsize=11,
                    anno_xshift=0.01,
                    marker_size=(100,100),
                    titles=["Male","Female"],
                    title_kwargs={"size":20},
                    anno_kwargs={"rotation":0,"fontsize":12}, verbose=False, check=False)
```

**stdout:**
```
2025/12/25 23:57:13  -POS data type is already numeric. Skipping conversion...
2025/12/25 23:57:13  -CHR data type is already numeric. Skipping conversion...
2025/12/25 23:57:24  -POS data type is already numeric. Skipping conversion...
2025/12/25 23:57:24  -CHR data type is already numeric. Skipping conversion...
```

```
(<Figure size 3200x2400 with 6 Axes>,
 <gwaslab.info.g_Log.Log at 0x7fd21de3fb30>)
```

![Output image](images/notebooks/visualization_stacked_mqq_img_6.png)

## Stacked regional plots with chromatin status

```python
gl.plot_stacked_mqq(objects=[gl1,gl2],
                    vcfs=[gl.get_path("1kg_eas_hg19"), gl.get_path("1kg_eas_hg19")],
                    region=(19,46214297 - 100000, 46214297 + 100000),
                    region_ref =["rs35560038","rs1055220","rs4802274"],
                    build="19",
                    anno="SNPID",
                    fontsize=15,
                    mode="r",
                    region_chromatin_files=["../0_sample_data/E098_15_coreMarks_mnemonics.bed.gz"],
                    region_chromatin_labels=["E098_15"],
                    cbar_fontsize=11,
                    anno_xshift=0.01,
                    marker_size=(100,100),
                    titles=["Male","Female"],
                    title_kwargs={"size":20},
                    anno_kwargs={"rotation":0,"fontsize":10}, 
                    verbose=False, check=False)
```

**stdout:**
```
2025/12/25 23:57:30  -POS data type is already numeric. Skipping conversion...
2025/12/25 23:57:30  -CHR data type is already numeric. Skipping conversion...
2025/12/25 23:57:40  -POS data type is already numeric. Skipping conversion...
2025/12/25 23:57:40  -CHR data type is already numeric. Skipping conversion...
```

```
(<Figure size 3200x3200 with 7 Axes>,
 <gwaslab.info.g_Log.Log at 0x7fd21de3fb30>)
```

![Output image](images/notebooks/visualization_stacked_mqq_img_7.png)

### turn off legend marker

```python
gl.plot_stacked_mqq(objects=[gl1,gl2],
                    vcfs=[gl.get_path("1kg_eas_hg19"), gl.get_path("1kg_eas_hg19")],
                    region_legend_marker=False,
                    region=(19,46214297 - 100000, 46214297 + 100000),
                    region_ref =["rs35560038","rs1055220","rs4802274"],
                    build="19",
                    anno="SNPID",
                    fontsize=15,
                    mode="r",
                    region_chromatin_files=["../0_sample_data/E098_15_coreMarks_mnemonics.bed.gz"],
                    region_chromatin_labels=["E098_15"],
                    cbar_fontsize=11,
                    marker_size=(100,100),
                    titles=["Male","Female"],
                    title_kwargs={"size":20},
                    anno_kwargs={"rotation":0,"fontsize":10}, 
                    verbose=False, check=False)
```

**stdout:**
```
2025/12/25 23:57:49  -POS data type is already numeric. Skipping conversion...
2025/12/25 23:57:49  -CHR data type is already numeric. Skipping conversion...
2025/12/25 23:57:59  -POS data type is already numeric. Skipping conversion...
2025/12/25 23:57:59  -CHR data type is already numeric. Skipping conversion...
```

```
(<Figure size 3200x3200 with 7 Axes>,
 <gwaslab.info.g_Log.Log at 0x7fd21de3fb30>)
```

![Output image](images/notebooks/visualization_stacked_mqq_img_8.png)
