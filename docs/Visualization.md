# Manhattan plot and QQ plot : plot_mqq()

GWASLab provides a customizable plotting function for Manhattan and Q-Q plots.

```
.plot_mqq()
```

## Examples

!!! example "Quick Manhattan and Q-Q plot without any options"
    ```python
    mysumstats.plot_mqq()
    ```
    
    <img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/196595098-23ff14eb-5579-4177-8d20-54816a410f48.png">

## Options

- [Using P or MLOG10P](#use-mlog10p-for-extreme-p-values)
- [Adjusting x axis](#x-axis-physical-position-or-rank)
- [Adjusting y axis](#y-axis-skip-low-and-shrink-high)
- [Changing layout](#manhattan-and-qq-plot-layout)
- [Annotation](#annotation)
- [Adding lines](#lines)
- [Highlight loci and Pinpoint variants](#highlight-specified-loci) 
- [Colors and fonts](#colors-and-fontsizes)
- [MAF-stratified QQ plot](#maf-stratified-qq-plot).
- [Changing titles](#titles)
- [Saving figures](#save-plots)

By setting the options, you can create highly customized Manhattan plots and Q-Q plots.

!!! example "A customized Manhattan and QQ plot"
    ```python
    mysumstats.plot_mqq(
                      mode="qqm",
                      cut=14,
                      skip=3, 
                      anno_set=["rs12509595","rs7989823"] ,
                      pinpoint=["rs7989823"], 
                      highlight=["rs12509595","19:15040733:T:C"],
                      highlight_windowkb =1000,
                      stratified=True,
                      marker_size=(5,10),
                      figargs={"figsize":(15,5),"dpi":300})
    ```
    
    <img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/196594621-840217aa-117d-49ac-ab58-15a5fa6675b1.png">

### Manhattan and QQ plot layout

- `mode` : determine the layout of manhattan plot and qq plot.
    - `mqq`: left manhatan, right QQ
    - `qqm`: left QQ , right manhatan
    - `"m"`: only manhattan plot
    - `"qq"`: only qq plot
    - `mqqratio`: width ratio

!!! info "Layout"
    <img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/196593277-0908d49e-40aa-4fe3-b214-d774ab4d0382.png">

----------------------------------------------------------------

### Use MLOG10P for extreme P values

- `scaled` : `boolean`. By default, GWASLab uses P values for mqq plot. But you can set `scaled=Ture` to use MLOG10P to plot.

!!! note "Variant with extreme P values"
    To plot the variant with extreme P values (P < 1e-300), you can use `scaled=False` to create the plot with MLOG10P instead of raw P values. To calculate MLOG10P for extreme P values from BETA/SE or Z scores, you can use `mysumstats.fill_data(to_fill=["MLOG10P"], extreme=True)`. For details, please refer to the "Extreme P values" section in [https://cloufield.github.io/gwaslab/Conversion/](https://cloufield.github.io/gwaslab/Conversion/).

### X axis: Physical position or rank

- `use_rank` : `boolean`. If True, GWASLab will use position rank instead of the physical base-pair positions for x aixs.

!!! note
    If using rank, there will be no gap in the plot. If using base-pair positions, certain regions of the chromosome might be reflected in the plot like the heterochromatin.

### Y axis: Skip "low" and shrink "high"

- `skip` : sometimes it is not necessary to plot all variants, we can skip the variants with low -log10(P) values. For example, we can exclude varints with -log10(P) lower than 3 from the plot by specifying `skip=3`
- `cut` : loci with extremly large -log10(P) value are very likely to dwarf other significant loci, so we want to scale down the -log10(P) for variants above a certain threshold. 
- `cutfactor`:  shrinkage factor, default is 10.
- `cut_line_color`: the color of the line above which y axis is rescaled.
- `sig_level=5e-8` : genome-wide significance threshold.

!!! info "Auxiliary lines"
    <img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/196593118-b7358559-71ef-49ae-b8b0-15b34031c142.png">

!!! note
    lambda GC calculation for QQ plot will not be affected by skip and cut. The calculation is conducted using all variants in the original dataset.

----------------------------------------------------------------

### Annotation

- `anno`: `boolean` or `string` or `"GENENAME"`   
    - `boolean`: `anno = True`,  the variants to annotate will be selected automatically using a sliding window with `windowsize=500`kb. chr:pos
    - `string`: the column name used for annotation
    - `"GENENAME"` : automatically annotate nrearest gene names, using pyensembl. (remember to specify `build`, default is `build="19"`)
- `repel_force` : when the annotation overlaps with other, try increasing the repel_force to increase the padding between annotations. (deault is 0.01)
- `anno_set `: if you want to annotate only a few specific variants, you can simply provide a list of SNPIDs. 
- `anno_alias` : snpid:text dictionary for customized annotation

!!! info "Repel force"
    
    <img width="300" alt="image" src="https://user-images.githubusercontent.com/40289485/196592902-36c52102-09f2-4c58-894b-a7177e87bd09.png">


!!! example "Skip variants with -log10P<3 and annotate the lead variants with chr:pos"

    ```python
    mysumstats.plot_mqq(skip=3,anno=True)
    ```

    <img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/196591270-592ce820-c541-4b84-a766-0b58bc6423ee.png">

!!! example "Skip variants with -log10P<3 and annotate the lead variants with GENENAME"
    ```python
    mysumstats.plot_mqq(skip=3,anno="GENENAME",build="19")
    ```

    <img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/196591371-262d31d5-9640-474f-af0d-d6c511c77280.png">

!!! example "Skip variants with -log10P<3 and annotate the variants in `anno_set`"
    ```python
    mysumstats.plot_mqq(skip=3, anno_set=["rs12509595","19:15040733:T:C"])
    ```

    <img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/196591966-c9618c45-456b-4eb8-991b-66420f847a97.png">

!!! example "Skip variants with -log10P<3 and annotate the variants in `anno_set` with alias in `anno_alias`" 
    ```python
    mysumstats.plot_mqq(skip=3, anno_set=["rs12509595","19:15040733:T:C"], anno_alias={"rs12509595":"anything you want here"})
    ```
    
    <img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/196592136-f2cbc488-e02f-409b-b0ac-2ab16f7b1fd4.png">

### Annotation style

Added since 3.3.23

GWASLab now support 3 types of annotation styles:

- `expand`
- `right`
- `tight`

!!! example "`anno_style="expand"`"
    
    <img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/215389210-f8b76b07-0b41-4932-9ba2-15f4a5508855.png">

!!! example "`anno_style="right"`"

    <img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/215389238-735c15f0-ea46-442a-b624-4e903c123b55.png">

!!! example "`anno_style="tight"`"

    <img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/215481722-22c8fdcc-b016-4277-ab5f-eb2924f61e01.png">


### Adjust arm position

- `anno_d`:`dict`,key is the number of arm starting form 0, value is the direction you want the arm to shift towards . For example, `anno_d = {4:"r"}` means shift the 4th arm to the right 
- `arm_offset`: `float` shift distance in points
- `arm_scale`: `float` factors to adjust the height for all arms
- `arm_scale_d`: `dict` factors to adjust the height for specific arms. key is the number of arm startinf form 0, value is the factor which will be multiplied to arm height.

!!! example "Adjust the direction the first to left and the thrd to right"
    ```python
    mysumstats.plot_mqq(skip=2,anno=True)
    ```
    
    <img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/197342763-ffd4b3c1-d57a-4351-8f42-fb91ae282d32.png">
    
    ```python
    mysumstats.plot_mqq(skip=2,anno=True,          
                        anno_d={1:"l",3:"r"},
                        arm_offset=50)
    ```
    
    <img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/197344117-2dc1261f-7784-48b2-9015-bdb0e73fce02.png">

!!! example "Adjust the length of arm"
    ```python
    mysumstats.plot_mqq(skip=2,anno=True,arm_scale=1.5)
    ```

    <img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/197345571-6a6ccb4e-6837-475d-ae3b-3c59fb447c56.png">

!!! example "Adjust the length of arm for each variant"
    ```python
    mysumstats.plot_mqq(skip=2,anno=True,arm_scale_d={1:1.5,2:1.2,3:1.1})
    ```
    
    <img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/197345530-08f4c9f6-d804-4e7e-865b-0b219179f6a9.png">

----------------------------------------------------------------

### Highlight specified loci

Highlight specified loci (color all variants in a region by specifying variants and the length of flanking regions). 

- `highlight ` : `list` specify the variants of loci for highlighting.
- `highlight_color`: `string` specify the color ussed for highlighting.
- `highlight_windowkb` : `int` specify the span of highlighted region (deault: `highlight_windowkb = 500` kp)

### Pinpoint specified variants

Pinpoint certain variants in the manhattan plot.

- `pinpoint` : a list of SNPIDs
- `pinpoint_color` : color for pinpoint

!!! example "Highlight loci and pinpoint variants"

    ```python
    mysumstats.plot_mqq(skip=3,anno="GENENAME",build="19",
                       highlight=["rs12509595","rs7989823"],
                       pinpoint=["rs671","19:15040733:T:C"])
    ```

    <img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/196593330-1794223c-79fd-40f9-942d-8acbbe00827b.png">

----------------------------------------------------------------

### Lines
- `sig_line`: boolean. if True, plot the significant threshold line (default: True)
- `sig_level`: float. The significance threshold (default:5e-8)
- `sig_level_lead`: float. The significance threshold for extracting lead variants to annotate (default:5e-8)
- `sig_line_color`:color. Significant threshold line color (default:"grey")
- `suggestive_sig_line`:boolean. if True, plot the suggestive line (default:False),
- `suggestive_sig_level`:float. The suggestive level (default:5e-6),
- `suggestive_sig_line_color`:color. Suggestive level line color (default:"grey"),
- `additional_line`: list of threshold. (default:None)
- `additional_line_color`: list of colors. (default:None)
- `cut_line_color`: color. the color for cut line (default: "#ebebeb")

!!! example "Plot lines"
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
                    additional_line_color=["yellow","green"])
    ```
    <img width="600" alt="image" src="https://github.com/Cloufield/gwaslab/assets/40289485/e0a8c92b-a3b5-430e-b11c-114f970ca34a">




### MAF-stratified QQ plot
- `stratified`: `boolean` if True, plot MAF straitified QQ plot. Require EAF in sumstats.
- `maf_bins`: `list` maf bins for straitification. (default: `maf_bins=[(0, 0.01), (0.01, 0.05), (0.05, 0.25),(0.25,0.5)]`)
- `maf_bin_colors`: `list` colors used for each bin. (default: `maf_bin_colors = ["#f0ad4e","#5cb85c", "#5bc0de","#000042"]`)

!!! example "MAF-stratified Q-Q plot"
    <img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/196593456-93539e3e-e8c0-476f-ae4e-e5843f970ba6.png">

----------------------------------------------------------------

### Colors and Fontsizes

```python
mysumstats.plot_mqq(
          colors=["#597FBD","#74BAD3"],
          cut_line_color="#ebebeb",
          sig_line_color="grey",
          highlight_color="#CB132D",
          pinpoint_color ="red",
          maf_bin_colors = ["#f0ad4e","#5cb85c", "#5bc0de","#000042"],
          fontsize = 10,
          anno_fontsize = 10,
          title_fontsize = 13,
          marker_size=(5,25)
)
```

- `colors` : `list`, a list of color for Manhattan plot. (default: `["#597FBD","#74BAD3"]`)
- `cut_line_color` : `string`, color for cut line.  (default: `"#ebebeb"`)
- `sig_line_color` : `string`, color for significance threshold line. (default: `"grey"`)
- `highlight_color` : `string`, color for highlighting. (default: `"#CB132D"`)
- `pinpoint_color` : `string`, color for pinpointing. (default: `"red"`)
- `maf_bin_colors` : `list`, a list of color formaf stratified Q-Q plot. (default: `"["#f0ad4e","#5cb85c", "#5bc0de","#000042"]"`)
- `fontsize = 10`: fontsize for ticklabels.
- `title_fontsize = 13`, fontsize for title.
- `anno_fontsize = 10`, fontsize for annotation.

!!! example
    ```python
    mysumstats.plot_mqq(skip=2,
                        cut=20,
                        colors=sns.color_palette("Set3"),
                        sig_line_color="red",
                        fontsize = 8)
    ```
    <img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/217693629-98261192-c521-4466-a9eb-6e9a4d4bd2cb.png">


### Titles

```python
mysumstats.plot_mqq(
          title =None,
          mtitle=None,
          qtitle=None,
          title_pad=1.08
        )
```

- title : `string` , title for the whole plot
- mtitle : `string` , title for the Manhattan plot
- qtitle : `string` , title for the Q-Q plot
- title_pad : `float` , padding for title

<img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/197341807-e5fb02da-64ba-4320-b001-bb9f5760fdc5.png">

----------------------------------------------------------------

### Figure settings

```
figargs= {"figsize":(15,5),"dpi":100}
```

- `figargs` : `dict`, key-values pairs that are passed to matplotlib `plt.subplots()`
    - `figsize` : figure size
    - `dpi` : dots per inch. For pulications, dpi>=300 is on of the common criteria.


### Save plots

```python
mysumstats.plot_mqq(save="mymqqplots.png",saveargs={"dpi":400,"facecolor":"white"})
```
Two options for saving plots in `.plot_mqq`

- `save`:`string`,path for saved plot
- `saveargs`: `dict`, other parameters passed to matplotlib `savefig` function.

!!! example

    - save as png: `mysumstats.plot_mqq(save="mymqqplots.png",saveargs={"dpi":300})`
    - save as PDF: `mysumstats.plot_mqq(save="mymqqplots.pdf",saveargs={"dpi":300})`
