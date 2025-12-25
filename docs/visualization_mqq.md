# Manhattan and Q-Q plot

## Load gwaslab 

```python
import gwaslab as gl
```

## Load data into Sumstats Object

```python
mysumstats = gl.Sumstats("../0_sample_data/t2d_bbj.txt.gz",
             snpid="SNP",
             chrom="CHR",
             pos="POS",
             ea="ALT",
             nea="REF",
             neaf="Frq",
             p="P",
             build="19",             
             verbose=False)
mysumstats.fix_chr(verbose=False)
```

```
<gwaslab.g_Sumstats.Sumstats at 0x7f605a338b00>
```

## Create Manhattan plot and QQ plot

```python
mysumstats.plot_mqq(skip=2)
```

**stdout:**
```
2025/12/25 23:38:18 Configured plot style for plot_mqq:mqq
2025/12/25 23:38:18 Starting Manhattan-QQ plot creation (Version v4.0.0)
2025/12/25 23:38:18  -Genomic coordinates are based on GRCh37/hg19...
2025/12/25 23:38:18  - Genomic coordinates version: 19 ...
2025/12/25 23:38:18  - Genome-wide significance level to plot is set to 5e-08 ...
2025/12/25 23:38:18  - Input sumstats contains 12557761 variants...
2025/12/25 23:38:18  - Manhattan-QQ plot layout mode selected: mqq
2025/12/25 23:38:19  -POS data type is already numeric. Skipping conversion...
2025/12/25 23:38:20  -CHR data type is already numeric. Skipping conversion...
2025/12/25 23:38:20 Finished loading specified columns from the statistics
2025/12/25 23:38:20 Start data conversion and sanity check:
2025/12/25 23:38:25  -Sumstats P values are being converted to -log10(P)...
2025/12/25 23:38:26  -Converting data above cut line...
2025/12/25 23:38:26  -Maximum -log10(P) value is 167.58838029403677 .
2025/12/25 23:38:26 Finished data conversion and sanity check.
2025/12/25 23:38:26 Start to create Manhattan-QQ plot with 332882 variants...
2025/12/25 23:38:26  -Creating background plot...
2025/12/25 23:38:27 Finished creating Manhattan-QQ plot successfully
2025/12/25 23:38:27 Start to extract variants for annotation...
2025/12/25 23:38:27  -Found 89 significant variants with a sliding window size of 500 kb...
2025/12/25 23:38:27 Finished extracting variants for annotation...
2025/12/25 23:38:27 Start to create QQ plot with 332882 variants:
2025/12/25 23:38:27  -Plotting all variants...
2025/12/25 23:38:27  -Expected range of P: (0,1.0)
2025/12/25 23:38:28  -Lambda GC (MLOG10P mode) at 0.5 is   1.21283
2025/12/25 23:38:28 Finished creating QQ plot successfully!
2025/12/25 23:38:28 Finished creating plot successfully
```

![Output image](images/notebooks/visualization_mqq_img_0.png)

## gl.plot_mqq() Options

### Layout mode

If plotting all variants, it may take several minutes. You can use `skip` to skip variants with low MLOG10P in the plot. 

Note: use verbose=False to stop printing log and use check=False to skip sanity check for mqq plots

4 patterns of layout: 

- `mode= "mqq"` (default)
- `mode= "qqm"`
- `mode= "qq"` 
- `mode= "m"`

```python
mysumstats.plot_mqq(skip=3,cut=20,mode="qqm",check=False, verbose=False)
```

**stdout:**
```
2025/12/25 23:38:32  -POS data type is already numeric. Skipping conversion...
2025/12/25 23:38:34  -CHR data type is already numeric. Skipping conversion...
```

![Output image](images/notebooks/visualization_mqq_img_1.png)

```python
mysumstats.plot_mqq(skip=3,cut=20,mode="qqm",check=False, verbose=False)
```

**stdout:**
```
2025/12/25 23:38:40  -POS data type is already numeric. Skipping conversion...
2025/12/25 23:38:42  -CHR data type is already numeric. Skipping conversion...
```

![Output image](images/notebooks/visualization_mqq_img_2.png)

```python
mysumstats.plot_mqq(skip=3,cut=20,mode="m",check=False,verbose=False)
```

**stdout:**
```
2025/12/25 23:38:48  -POS data type is already numeric. Skipping conversion...
2025/12/25 23:38:49  -CHR data type is already numeric. Skipping conversion...
```

![Output image](images/notebooks/visualization_mqq_img_3.png)

```python
mysumstats.plot_mqq(mode="qq", fig_args= {"figsize":(1,1)},check=False,verbose=False)
```

![Output image](images/notebooks/visualization_mqq_img_4.png)

### Y axis

#### skip

- `skip` : skip the variants with low -log10(P) values for plotting


```python
mysumstats.plot_mqq(skip=5,mode= "m", check=False ,verbose=False)
```

**stdout:**
```
2025/12/25 23:39:30  -POS data type is already numeric. Skipping conversion...
2025/12/25 23:39:32  -CHR data type is already numeric. Skipping conversion...
```

![Output image](images/notebooks/visualization_mqq_img_5.png)

#### cut

 `cut` : scale down the -log10(P) for variants above a certain threshold

```python
mysumstats.plot_mqq(cut=20 ,mode= "m",check=False,verbose=False)
```

**stdout:**
```
2025/12/25 23:39:36  -POS data type is already numeric. Skipping conversion...
2025/12/25 23:39:38  -CHR data type is already numeric. Skipping conversion...
```

![Output image](images/notebooks/visualization_mqq_img_6.png)

Make the Y axis jagged to indicate that it has been rescale. 

#### jagged

```python
mysumstats.plot_mqq(skip=3,cut=20, jagged=True,mode= "m",check=False,verbose=False)
```

**stdout:**
```
2025/12/25 23:41:06  -POS data type is already numeric. Skipping conversion...
2025/12/25 23:41:07  -CHR data type is already numeric. Skipping conversion...
```

![Output image](images/notebooks/visualization_mqq_img_7.png)

### X axis

#### use_rank

- `use_rank`: if True, GWASLab will use position rank instead of the physical base-pair positions for x aixs.

There will be no gap if `use_rank = True`

```python
mysumstats.plot_mqq(skip=3,cut=20,use_rank=True,mode= "m",check=False,verbose=False)
```

**stdout:**
```
2025/12/25 23:41:12  -POS data type is already numeric. Skipping conversion...
2025/12/25 23:41:13  -CHR data type is already numeric. Skipping conversion...
```

![Output image](images/notebooks/visualization_mqq_img_8.png)

#### xtight

`xtight=True` can be used to remove the padding 

```python
mysumstats.plot_mqq(xtight=True, skip=3,cut=20,mode= "m",check=False,verbose=False)
```

**stdout:**
```
2025/12/25 23:41:18  -POS data type is already numeric. Skipping conversion...
2025/12/25 23:41:19  -CHR data type is already numeric. Skipping conversion...
```

![Output image](images/notebooks/visualization_mqq_img_9.png)

#### chrpad

`chrpad`: adjust space between each chromosome by `max(POS) * chrpad`

```python
mysumstats.plot_mqq(chrpad=0.2,xtight=True, skip=3,cut=20,mode= "m",check=False,verbose=False)
```

**stdout:**
```
2025/12/25 23:41:23  -POS data type is already numeric. Skipping conversion...
2025/12/25 23:41:25  -CHR data type is already numeric. Skipping conversion...
```

![Output image](images/notebooks/visualization_mqq_img_10.png)

### Annotation

#### anno=True

`anno=True` : annoatate all lead variants with chr:pos

```python
mysumstats.plot_mqq(skip=3,cut=20,anno=True,mode= "m",check=False,verbose=False)
```

**stdout:**
```
2025/12/25 23:41:29  -POS data type is already numeric. Skipping conversion...
2025/12/25 23:41:31  -CHR data type is already numeric. Skipping conversion...
```

![Output image](images/notebooks/visualization_mqq_img_11.png)

Since there are a large number of novel loci, if we annotate all loci, it will be too messy. Let's only annotate the loci with P<1e-20 by specifying `sig_level_lead=1e-20`.

```python
mysumstats.plot_mqq(skip=3,cut=20,anno=True, anno_sig_level=1e-40,mode= "m",check=False,verbose=True)
```

**stdout:**
```
2025/12/25 23:41:36 Configured plot style for plot_mqq:m
2025/12/25 23:41:36 Starting Manhattan plot creation (Version v4.0.0)
2025/12/25 23:41:36  -Genomic coordinates are based on GRCh37/hg19...
2025/12/25 23:41:36  - Genomic coordinates version: 19 ...
2025/12/25 23:41:36  - Genome-wide significance level to plot is set to 5e-08 ...
2025/12/25 23:41:36  - Input sumstats contains 12557761 variants...
2025/12/25 23:41:36  - Manhattan plot layout mode selected: m
2025/12/25 23:41:36  -POS data type is already numeric. Skipping conversion...
2025/12/25 23:41:38  -CHR data type is already numeric. Skipping conversion...
2025/12/25 23:41:38 Finished loading specified columns from the statistics
2025/12/25 23:41:38 Start data conversion and sanity check:
2025/12/25 23:41:38  -Sanity check will be skipped.
2025/12/25 23:41:39  -Sumstats P values are being converted to -log10(P)...
2025/12/25 23:41:41  -Converting data above cut line...
2025/12/25 23:41:41  -Maximum -log10(P) value is 167.58838029403677 .
2025/12/25 23:41:41  -Minus log10(P) values above 20 will be shrunk with a shrinkage factor of 10...
2025/12/25 23:41:41 Finished data conversion and sanity check.
2025/12/25 23:41:41 Start to create Manhattan plot with 91234 variants...
2025/12/25 23:41:41  -Creating background plot...
2025/12/25 23:41:41 Finished creating Manhattan plot successfully
2025/12/25 23:41:41 Start to extract variants for annotation...
2025/12/25 23:41:41  -Found 7 significant variants with a sliding window size of 500 kb...
2025/12/25 23:41:41 Finished extracting variants for annotation...
2025/12/25 23:41:41  -Annotating using column CHR:POS...
2025/12/25 23:41:41  -Adjusting text positions with repel_force=0.03...
2025/12/25 23:41:41 Finished creating plot successfully
```

![Output image](images/notebooks/visualization_mqq_img_12.png)

#### anno="GENENAME"

- `anno="GENENAME"` : automatically annoatate the nearest gene name

Note: remerber to set `build=19` or `build=38` when loading or plotting.

```python
mysumstats.plot_mqq(skip=3,cut=20,anno="GENENAME", sig_level_lead=1e-20,mode= "m",check=False,verbose=False)
```

**stdout:**
```
2025/12/25 23:41:42  -POS data type is already numeric. Skipping conversion...
2025/12/25 23:41:44  -CHR data type is already numeric. Skipping conversion...
```

![Output image](images/notebooks/visualization_mqq_img_13.png)

You can specify `anno_gtf_path` to use you own GTF file for GENENAME annotation

```python
# mysumstats.plot_mqq(anno_gtf_path="/home/yunye/.gwaslab/Homo_sapiens.GRCh37.87.chr.gtf.gz",anno="GENENAME")
```

#### anno_d

We can use `anno_d` to slightly adjust the arrows.

`anno_d` accepts a dictionary of `index of annotation`: `left`/`right`

For example, `1:"left"` means to adjust towards left. 

```python
mysumstats.plot_mqq(anno_d={1:"left",2:"right"}, skip=3,cut=20,anno=True, sig_level_lead=1e-20,mode= "m",check=False,verbose=False)
```

**stdout:**
```
2025/12/25 23:41:53  -POS data type is already numeric. Skipping conversion...
2025/12/25 23:41:54  -CHR data type is already numeric. Skipping conversion...
```

![Output image](images/notebooks/visualization_mqq_img_14.png)

#### anno_scale

We can also use `arm_scale` to adjust where to put the annotation texts.

For example, `arm_scale=1.2` means the default length will be multiplied by a factor of 1.2.

```python
mysumstats.plot_mqq(arm_scale=1.2, skip=3,cut=20,anno=True, sig_level_lead=1e-20,mode= "m",check=False,verbose=False)
```

**stdout:**
```
2025/12/25 23:41:59  -POS data type is already numeric. Skipping conversion...
2025/12/25 23:42:01  -CHR data type is already numeric. Skipping conversion...
```

![Output image](images/notebooks/visualization_mqq_img_15.png)

`arm_scale_d` accepts a dictionary of `index of annotation`: `arm_scale`

For example, `1:1.2` means to adjust the arm of the second by a factor of 1.2. 

```python
mysumstats.plot_mqq(arm_scale_d={0:0.8,
                                 1:0.7,
                                 2:0.6,
                                 3:0.8}, 
                    skip=3,cut=20,anno=True, sig_level_lead=1e-20,mode= "m",check=False,verbose=False)
```

**stdout:**
```
2025/12/25 23:42:07  -POS data type is already numeric. Skipping conversion...
2025/12/25 23:42:08  -CHR data type is already numeric. Skipping conversion...
```

![Output image](images/notebooks/visualization_mqq_img_16.png)

#### anno_style

GWASLab provides three types of different annotation styles

`anno_style="right"`, `anno_style="expand"`, and `anno_style="tight"`

```python
mysumstats.plot_mqq(skip=3,cut=20,anno=True,anno_style="expand", sig_level_lead=1e-20,mode= "m",check=False,verbose=False)
```

**stdout:**
```
2025/12/25 23:42:14  -POS data type is already numeric. Skipping conversion...
2025/12/25 23:42:15  -CHR data type is already numeric. Skipping conversion...
```

![Output image](images/notebooks/visualization_mqq_img_17.png)

```python
mysumstats.plot_mqq(skip=3,cut=20,anno=True,anno_style="tight", sig_level_lead=1e-20,mode= "m", check=False,verbose=False)
```

**stdout:**
```
2025/12/25 23:42:21  -POS data type is already numeric. Skipping conversion...
2025/12/25 23:42:22  -CHR data type is already numeric. Skipping conversion...
```

![Output image](images/notebooks/visualization_mqq_img_18.png)

#### anno_set

If we want to annotate only a subset of variants, we can pass a list of variant IDs to `anno_set`. 

Let's check all lead variants and select only two to annotate.

```python
mysumstats.get_lead(verbose=False).sort_values(by="P")
```

```
| SNPID | CHR | POS | EA NEA | STATUS | EAF | \ |
| --- | --- | --- | --- | --- | --- | --- |
| 7947442 | 11:2858546_C_T | 11 | 2858546 | C | T | 1995999 |
| 6903480 | 9:22132729_A_G | 9 | 22132729 | G | A | 1995999 |
| 4664126 | 6:20688121_T_A | 6 | 20688121 | T | A | 1995999 |
| 5965364 | 7:127253550_C_T | 7 | 127253550 | C | T | 1995999 |
| 12554976 | X:152908887_G_A | 23 | 152908887 | G | A | 1995999 |
| ... | ... | ... | ... | ... | ... | ... |
| 12272930 | 2.616000e-08 |  |  |  |  |  |
| 4599424 | 2.849000e-08 |  |  |  |  |  |
| 10287151 | 3.134000e-08 |  |  |  |  |  |
| 534095 | 3.289000e-08 |  |  |  |  |  |
| 10874150 | 4.159000e-08 |  |  |  |  |  |

*[89 rows x 8 columns]*
```

This time, let's annotate 1:154309595_TA_T and 2:27734972_G_A with its nearest gene names!


`anno_set` : the set of variants you want to annotate

```python
mysumstats.plot_mqq(skip=3,cut=20,anno="GENENAME",anno_set=["1:154309595_TA_T","2:27734972_G_A"],mode= "m",check=False,verbose=False)
```

**stdout:**
```
2025/12/25 23:42:55  -POS data type is already numeric. Skipping conversion...
2025/12/25 23:42:56  -CHR data type is already numeric. Skipping conversion...
```

![Output image](images/notebooks/visualization_mqq_img_19.png)

#### anno_alias

anno_alias accepts a dictionary of `SNPID`:`string`. You can use this to customized the text for annotation.

```python
mysumstats.plot_mqq(skip=3,
                    cut=20,
                    anno=True,
                    anno_set=["1:154309595_TA_T","2:27734972_G_A"],
                    anno_alias={"1:154309595_TA_T":"anything you want here"},
                    mode= "m",
                    verbose=False)
```

**stdout:**
```
2025/12/25 23:43:03  -POS data type is already numeric. Skipping conversion...
2025/12/25 23:43:04  -CHR data type is already numeric. Skipping conversion...
```

![Output image](images/notebooks/visualization_mqq_img_20.png)

### Highlight loci & Pinpoint variants (single group)

- `highlight`: a variant list of loci you want to highlight
- `pinpoint`: a variant list of variants you want to pinpoint

```python
mysumstats.plot_mqq(skip=3,cut=20,
                    anno=True,
                    anno_set=["2:27734972_G_A","2:27734972_G_A","7:127253550_C_T", "19:46166604_C_T"],
                    highlight=["19:46166604_C_T","1:154309595_TA_T","7:127253550_C_T"],
                    highlight_windowkb=1000,
                    pinpoint=["2:27734972_G_A"],
                    mode= "m",
                    check=False,verbose=False)
```

**stdout:**
```
2025/12/25 23:43:11  -POS data type is already numeric. Skipping conversion...
2025/12/25 23:43:13  -CHR data type is already numeric. Skipping conversion...
```

![Output image](images/notebooks/visualization_mqq_img_21.png)

### Highlight loci & Pinpoint variants (multi-group)

Instead of a list, you can provide a list of lists. Each member list is then a group.

```python
mysumstats.plot_mqq(skip=3,cut=20,
                    highlight=[
                               ["19:46166604_C_T","1:154309595_TA_T"],
                               ["X:57170781_A_AT","7:127253550_C_T"]
                              ],
                    highlight_windowkb=1000,
                    highlight_color=["yellow","red"],
                    mode= "m",
                    check=False,verbose=False)
```

**stdout:**
```
2025/12/25 23:43:20  -POS data type is already numeric. Skipping conversion...
2025/12/25 23:43:22  -CHR data type is already numeric. Skipping conversion...
```

![Output image](images/notebooks/visualization_mqq_img_22.png)

### MAF-stratified QQ plot

```python
mysumstats.plot_mqq(skip=3,cut=20, mode="mqq",stratified=True,check=False,verbose=False)
```

**stdout:**
```
2025/12/25 23:43:29  -POS data type is already numeric. Skipping conversion...
2025/12/25 23:43:30  -CHR data type is already numeric. Skipping conversion...
2025/12/25 23:43:30  -EAF data type is already numeric. Skipping conversion...
```

![Output image](images/notebooks/visualization_mqq_img_23.png)

### Auxiliary lines

```python
mysumstats.plot_mqq(skip=3,
                build="19",
                anno="GENENAME",
                windowsizekb=1000000,
                cut=20,
                cut_line_color="purple",
                sig_level=5e-8,  
                sig_level_lead=1e-6, 
                sig_line_color="grey",
                suggestive_sig_line = True,
                suggestive_sig_level = 1e-6,
                suggestive_sig_line_color="blue",
                additional_line=[1e-40,1e-60],
                additional_line_color=["yellow","green"],
                mode= "m",check=False,verbose=False)
```

**stdout:**
```
2025/12/25 23:43:38  -POS data type is already numeric. Skipping conversion...
2025/12/25 23:43:40  -CHR data type is already numeric. Skipping conversion...
```

![Output image](images/notebooks/visualization_mqq_img_24.png)

### Font and marker size

- `fontsize`
- `anno_fontsize`
- `title_fontsize`
- `marker_size`

```python
mysumstats.plot_mqq(skip=3,
          anno="GENENAME",
          title= "My Manhattan Plot",
          windowsizekb=1000000,
          fontsize =13,
          anno_fontsize = 15,
          title_fontsize = 30,
          marker_size=(5,25),
          mode= "m", 
          verbose=False,
          check=False
)
```

**stdout:**
```
2025/12/25 23:43:47  -POS data type is already numeric. Skipping conversion...
2025/12/25 23:43:48  -CHR data type is already numeric. Skipping conversion...
```

![Output image](images/notebooks/visualization_mqq_img_25.png)

### Colors

- `colors`
- `cut_line_color`
- `sig_line_color`
- `highlight_color`
- `pinpoint_color`
- `maf_bin_colors`

```python
mysumstats.plot_mqq(skip=3,
          cut=20,
          stratified=True,
          highlight=["7:127253550_C_T"],
          pinpoint=["2:27734972_G_A"],
          colors=["orange","blue"],
          cut_line_color="yellow",
          sig_line_color="red",
          highlight_color="purple",
          pinpoint_color ="green",
          maf_bin_colors = ["#FFE2D1","#E1F0C4", "#6BAB90","#55917F"],
          check=False,verbose=False
)
```

**stdout:**
```
2025/12/25 23:43:55  -POS data type is already numeric. Skipping conversion...
2025/12/25 23:43:56  -CHR data type is already numeric. Skipping conversion...
2025/12/25 23:43:56  -EAF data type is already numeric. Skipping conversion...
```

![Output image](images/notebooks/visualization_mqq_img_26.png)

### Save plots

```python
mysumstats.plot_mqq(skip=3,
                    cut=20, 
                    mode="mqq",
                    stratified=True,
                    save="my_maf_stratified_mqq_plot.png",
                    save_kwargs={"dpi":300,"facecolor":"white"})
```

**stdout:**
```
2025/12/25 23:44:06 Configured plot style for plot_mqq:mqq
2025/12/25 23:44:06 Starting Manhattan-QQ plot creation (Version v4.0.0)
2025/12/25 23:44:06  -Genomic coordinates are based on GRCh37/hg19...
2025/12/25 23:44:06  - Genomic coordinates version: 19 ...
2025/12/25 23:44:06  - Genome-wide significance level to plot is set to 5e-08 ...
2025/12/25 23:44:06  - Input sumstats contains 12557761 variants...
2025/12/25 23:44:06  - Manhattan-QQ plot layout mode selected: mqq
2025/12/25 23:44:06  -POS data type is already numeric. Skipping conversion...
2025/12/25 23:44:08  -CHR data type is already numeric. Skipping conversion...
2025/12/25 23:44:08  -EAF data type is already numeric. Skipping conversion...
2025/12/25 23:44:08 Finished loading specified columns from the statistics
2025/12/25 23:44:08 Start data conversion and sanity check:
2025/12/25 23:44:13  -Removed 0 variants with nan in EAF column ...
2025/12/25 23:44:14  -Sumstats P values are being converted to -log10(P)...
2025/12/25 23:44:16  -Converting data above cut line...
2025/12/25 23:44:16  -Maximum -log10(P) value is 167.58838029403677 .
2025/12/25 23:44:16  -Minus log10(P) values above 20 will be shrunk with a shrinkage factor of 10...
2025/12/25 23:44:16 Finished data conversion and sanity check.
2025/12/25 23:44:16 Start to create Manhattan-QQ plot with 91234 variants...
2025/12/25 23:44:16  -Creating background plot...
2025/12/25 23:44:17 Finished creating Manhattan-QQ plot successfully
2025/12/25 23:44:17 Start to extract variants for annotation...
2025/12/25 23:44:17  -Found 89 significant variants with a sliding window size of 500 kb...
2025/12/25 23:44:17 Finished extracting variants for annotation...
2025/12/25 23:44:17 Start to create QQ plot with 91234 variants:
2025/12/25 23:44:17  -Plotting variants stratified by MAF...
2025/12/25 23:44:18  -Lambda GC (MLOG10P mode) at 0.5 is   1.21283
2025/12/25 23:44:18 Finished creating QQ plot successfully!
2025/12/25 23:44:18 Start to save figure...
2025/12/25 23:44:20  -Saved to my_maf_stratified_mqq_plot.png successfully! (png) (overwrite)
2025/12/25 23:44:20 Finished saving figure...
2025/12/25 23:44:20 Finished creating plot successfully
```

![Output image](images/notebooks/visualization_mqq_img_27.png)

### Fig object

plot_mqq will return a matplotlib figure object

```python
my_mqqplot = mysumstats.plot_mqq(skip=3,cut=20, mode="mqq",stratified=True,verbose=False)
```

**stdout:**
```
2025/12/25 23:44:41  -POS data type is already numeric. Skipping conversion...
2025/12/25 23:44:42  -CHR data type is already numeric. Skipping conversion...
2025/12/25 23:44:42  -EAF data type is already numeric. Skipping conversion...
```

```python
my_mqqplot
```

![Output image](images/notebooks/visualization_mqq_img_28.png)
