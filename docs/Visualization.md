# Manhattan plot and QQ plot 

GWASLab provides a customizable plotting function for Manhattan and Q-Q plots.

## .plot_mqq()

```
.plot_mqq()
```

## A simple example

!!! example "Quick Manhattan and Q-Q plot without any options"
    ```python
    mysumstats.plot_mqq()
    ```
    
    <img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/196595098-23ff14eb-5579-4177-8d20-54816a410f48.png">

See other examples [here](https://cloufield.github.io/gwaslab/visualization_mqq/).

## Options

- [Using P or MLOG10P](#22-use-mlog10p-for-extreme-p-values)
- [Adjusting x axis](#23-x-axis-physical-position-or-rank)
- [Adjusting y axis](#24-y-axis-skip-low-and-shrink-high)
- [Changing layout](#21-manhattan-and-qq-plot-layout)
- [Annotation](#25-annotation)
- [Adding lines](#210-lines)
- [Highlight loci and Pinpoint variants](#28-highlight-specified-loci) 
- [Colors and fonts](#212-colors-and-fontsizes)
- [MAF-stratified QQ plot](#211-maf-stratified-qq-plot).
- [Changing titles](#213-titles)
- [Saving figures](#215-saving-plots)

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

### Plot layout


|Option|DataType|Description|Default|
|-|-|-|-|
|`mode`|`mqq`,`qqm`,`qq`,`m`|Determine the layout of manhattan plot and qq plot. <br/> `mqq`: left manhatan, right QQ plot <br/>`qqm`: left QQ plot, right Manhattan plot <br/>`m`: only Manhattan plot<br/> `"qq"`: only qq plot|`mqq`|
|`mqqratio`|`float`|width ratio of Manhattan plot and QQ plot|`3`|

!!! info "Layout"
    <img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/196593277-0908d49e-40aa-4fe3-b214-d774ab4d0382.png">

----------------------------------------------------------------

### Use MLOG10P for extreme P values

|Option|DataType|Description|Default|
|-|-|-|-|
|`scaled`|`boolean`|By default, GWASLab uses P values for mqq plot. But you can set `scaled=Ture` to use MLOG10P to plot.|`False`|

!!! note "Variant with extreme P values"
    To plot the variant with extreme P values (P < 1e-300), you can use `scaled=False` to create the plot with MLOG10P instead of raw P values. To calculate MLOG10P for extreme P values from BETA/SE or Z scores, you can use `mysumstats.fill_data(to_fill=["MLOG10P"], extreme=True)`. For details, please refer to the "Extreme P values" section in [https://cloufield.github.io/gwaslab/Conversion/](https://cloufield.github.io/gwaslab/Conversion/).

### X axis: Physical position or rank

|Option|DataType|Description|Default|
|-|-|-|-|
|`use_rank`|`boolean`|If True, GWASLab will use position rank instead of the physical base-pair positions for x aixs.|`False`|

!!! note
    If using rank, there will be no gap in the plot. If using base-pair positions, certain regions of the chromosome might be reflected in the plot like the heterochromatin.

### Y axis: Skip "low" and shrink "high"

|Option|DataType|Description|Default|
|-|-|-|-|
|`skip`|`float`|Sometimes it is not necessary to plot all variants, we can skip the variants with low -log10(P) values for plotting. For example, we can omit varints with -log10(P) lower than 3 from the plot by specifying `skip=3`. Calculation of lambda GC won't be affected by this|`None`|
|`cut`|`float`|loci with extremly large -log10(P) value are very likely to dwarf other significant loci, so we want to scale down the -log10(P) for variants above a certain threshold. |`None`|
|`cutfactor`|`float`|shrinkage factor|`10`|
|`cut_line_color`|`float`|the color of the line above which y axis is rescaled.|`500`|
|`sig_level`|`float`|genome-wide significance threshold |`5e-8`|


!!! info "Auxiliary lines"
    <img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/196593118-b7358559-71ef-49ae-b8b0-15b34031c142.png">

!!! note
    lambda GC calculation for QQ plot will not be affected by skip and cut. The calculation is conducted using all variants in the original dataset.

----------------------------------------------------------------

### Annotation

|Option|DataType|Description|Default|
|-|-|-|-|
|`anno`|`boolean` or `string` or `"GENENAME"`|If `anno = True`, variants will anotated with chr:pos; or `string`, the column name used for annotation; or `"GENENAME"`, automatically annotate nrearest gene names, using pyensembl. (remember to specify `build`, default is `build="19"`) |`False`|
|`anno_set`|`list`|If you want to annotate only a few specific variants, you can simply provide a list of SNPIDs or rsIDs for annotation. If None, the variants to annotate will be selected automatically using a sliding window with `windowsize=500`kb.|`None`|
|`repel_force`|`float`|when the annotation overlaps with other, try increasing the repel_force to increase the padding between annotations.|`0.01`|
|`anno_alias`|`dict`|snpid:text dictionary for customized annotation|`None`|


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


### Adjust arm positions

|Option|DataType|Description|Default|
|-|-|-|-|
|`anno_d`|`dict`|key is the number of arm starting form 0, value is the direction you want the arm to shift towards . For example, `anno_d = {4:"r"}` means shift the 4th arm to the right |`None`|
|`arm_offset`|`float`|distance in points|`500`|
|`arm_scale`|`float`|factors to adjust the height for all arms|`1.0`|
|`arm_scale_d`|`dict`|factors to adjust the height for specific arms. key is the number of arm startinf form 0, value is the factor which will be multiplied to arm height.|`None`|

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

### Highlight loci

Highlight specified loci (color all variants in a region by specifying variants and the length of flanking regions). 

|Highlighting Option|DataType|Description|Default|
|-|-|-|-|
|`highlight`|`list`|a list of SNPID or rsID; these loci (all variants in the specified variants positions +/- `highlight_windowkb`) will be highlighted in `pinpoint_color`|`True`|
|`highlight_windowkb`|`int`|Specify the span of highlighted region in kbp|`500`|
|`highlight_color`|`list`|Color for highlighting loci|`"#CB132D"`|

### Pinpoint variants

Pinpoint certain variants in the Manhattan plot.

|Pinpointing Option|DataType|Description|Default|
|-|-|-|-|
|`pinpoint`|`list`|a list of SNPID or rsID; these variants will be highlighted in `pinpoint_color`|`True`|
|`pinpoint_color`|`list`|color for pinpointing variants|`"red"`|

!!! example "Highlight loci and pinpoint variants"

    ```python
    mysumstats.plot_mqq(skip=3,anno="GENENAME",build="19",
                       highlight=["rs12509595","rs7989823"],
                       pinpoint=["rs671","19:15040733:T:C"])
    ```

    <img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/196593330-1794223c-79fd-40f9-942d-8acbbe00827b.png">

----------------------------------------------------------------

### Lines

|Line Option|DataType|Description|Default|
|-|-|-|-|
|`sig_line`|`boolean`|If True, plot the significant threshold line|`True`|
|`sig_level`|`float`|The significance threshold|`5e-8`|
|`sig_level_lead`|`float`|The significance threshold for extracting lead variants to annotate|`5e-8`|
|`sig_line_color`|`string`|If True, plot the significant threshold line|`True`|
|`suggestive_sig_line`|`boolean`|If True, plot the suggestive threshold line|`True`|
|`suggestive_sig_level`|`float`|The suggestive threshold|`5e-6`|
|`suggestive_sig_line_color`|`string`|Suggestive level line color|`"grey"`|
|`additional_line`|`list`|list of P values used to plot additional lines|`None`|
|`additional_line_color`|`list`|list of colors for the additional lines|`None`|
|`cut_line_color`|`string`|If True, plot the significant threshold line|`"#ebebeb"`|

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

| QQ plot Option   | DataType  | Description                                                     | Default                                              |
|------------------|-----------|-----------------------------------------------------------------|------------------------------------------------------|
| `stratified`     | `boolean` | if True, plot MAF straitified QQ plot. Require EAF in sumstats. | `False`                                              |
| `maf_bins`       | `list`    | MAF bins for straitification.                                   | `[(0, 0.01), (0.01, 0.05), (0.05, 0.25),(0.25,0.5)]` |
| `maf_bin_colors` | `list`    | colors used for each MAF bin.                                   | `["#f0ad4e","#5cb85c", "#5bc0de","#000042"]`         |

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

Color-related options

| Color Option      | DataType | Description                                                                          | Default                                      |
|-------------------|----------|--------------------------------------------------------------------------------------|----------------------------------------------|
| `colors`          | `list`   | a list of colors for chromsomes in the Manhattan plot; it will be used repetitively. | `["#597FBD","#74BAD3"]`                      |
| `cut_line_color`  | `string` | color for the cut line.                                                              | `"#EBEBEB"`                                  |
| `sig_line_color`  | `string` | color for significance threshold line.                                               | `"grey"`                                     |
| `highlight_color` | `string` | color for highlighting loci                                                          | `"#CB132D"`                                  |
| `pinpoint_color`  | `string` | color for pinpointing variants                                                       | `"red"`                                      |
| `maf_bin_colors`  | `list`   | a list of colors for maf-stratified Q-Q plot.                                        | `["#f0ad4e","#5cb85c", "#5bc0de","#000042"]` |


Font-related options

| Font Option      | DataType | Description              | Default   |
|------------------|----------|--------------------------|-----------|
| `fontsize`       | `list`   | fontsize for ticklabels. | `9`       |
| `title_fontsize` | `13`     | fontsize for title.      | `13`      |
| `anno_fontsize`  | `10`     | fontsize for annotation. | `9`       |
| `font_family`    | `string` | font family              | `"Arial"` |

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

| Title Option | DataType | Description                  | Default |
|--------------|----------|------------------------------|---------|
| `title`      | `string` | title for the figure.        | ``      |
| `mtitle`     | `string` | title for the Manhattan plot | ``      |
| `qtitle`     | `string` | title for the Q-Q plot       | ``      |
| `title_pad`  | `float`  | padding for title            | `1.08`  |

<img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/197341807-e5fb02da-64ba-4320-b001-bb9f5760fdc5.png">

----------------------------------------------------------------

### Figure settings

```
figargs= {"figsize":(15,5),"dpi":100}
```

| Figure Option | DataType | Description                                                     | Default                        |
|---------------|----------|-----------------------------------------------------------------|--------------------------------|
| `figargs`     | `dict`   | key-values pairs that are passed to matplotlib `plt.subplots()` | `{"figsize":(15,5),"dpi":200}` |

Commonly used ones: 

- `figsize` : figure size
- `dpi` : dots per inch. For pulications, dpi>=300 is on of the common criteria.


### Saving plots

```python
mysumstats.plot_mqq(save="mymqqplots.png",save_args={"dpi":400,"facecolor":"white"})
```

Two options for saving plots in `.plot_mqq`

| Saving Option | DataType              | Description                                                                                            | Default                           |
|---------------|-----------------------|--------------------------------------------------------------------------------------------------------|-----------------------------------|
| `save`        | `string` or `boolean` | If `string`, the plot will be saved to the specified path; If `True`, it will be saved to default path | `True`                            |
| `save_args`   | `dict`                | other parameters passed to matplotlib `savefig` function.                                              | `{"dpi":300,"facecolor":"white"}` |

!!! example

    - save as png: `mysumstats.plot_mqq(save="mymqqplots.png",save_args={"dpi":300})`
    - save as PDF: `mysumstats.plot_mqq(save="mymqqplots.pdf",save_args={"dpi":300})`
