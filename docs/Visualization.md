## Manhattan plot and QQ plot : plot_mqq()

```
import gwaslab as gl

mydata = gl.Sumstats(....)

mydata.plot_mqq(
          mlog10p="MLOG10P",
          scaled=False,
          mode="mqq",
          mqqratio=3,
          windowsizekb=500,
          anno=None,
          anno_set=[],
          anno_alias={},
          anno_d={},
          arm_offset=50,
          arm_scale=1,
          cut=0,
          skip=0,
          cutfactor=10,
          cut_line_color="#ebebeb",  
          sig_level=5e-8,
          sig_line_color="grey",
          suggestive_sig_level=5e-6,
          highlight = [],
          highlight_color="#CB132D",
          highlight_windowkb = 500,
          pinpoint=[],
          pinpoint_color ="red",
          stratified=False,
          maf_bins=[(0, 0.01), (0.01, 0.05), (0.05, 0.25),(0.25,0.5)],
          maf_bin_colors = ["#f0ad4e","#5cb85c", "#5bc0de","#000042"],
          gc=True,
          title =None,
          mtitle=None,
          qtitle=None,
          figargs= {"figsize":(15,5),"dpi":100},
          fontsize = 10,
          colors=["#597FBD","#74BAD3"],
          marker_size=(5,25),
          use_rank=False,
          verbose=True,
          repel_force=0.03,
          build="19",
          title_pad=1.08, 
          save=None,
          saveargs={"dpi":400,"facecolor":"white"},
          ):
```

## Quick plot

```
mysumstatsysumstats.plot_mqq()
```

<img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/196595098-23ff14eb-5579-4177-8d20-54816a410f48.png">


## Customized plot

```
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



## Options

### Manhattan and QQ plot layout

`mode` : determine the layout of manhattan plot and qq plot.

`"mqq" `or `"qqm" `: side-by-side manhattan and QQ plt. 

- `mqq`: left manhatan, right QQ
- `qqm`: left QQ , right manhatan
- `"m"`: only manhattan plot
- `"qq"`: only qq plot
  
`mqqratio`: width ratio

<img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/196593277-0908d49e-40aa-4fe3-b214-d774ab4d0382.png">

----------------------------------------------------------------

### Skip "low" and shrink "high"

- `skip` : sometimes it is not necessary to plot all variants, we can skip the insignicant variants . For example, we can exclude varints with -log10p lower than 3 from the plot by specifying `skip=3`
- `cut` : loci with extremly large -log10(P) value are very likely to dwarf other significant loci , so we want to scale down the extrame loci from a certain threshold. 
- `cutfactor`:  shrinkage factor, default is 10 
- `cut_line_color`: the color of the line above which y axis is rescaled 

<img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/196593118-b7358559-71ef-49ae-b8b0-15b34031c142.png">

----------------------------------------------------------------

### Annotation
* `anno`: `boolean` or `string` or `"GENENAME"`   
    * `boolean`: `anno = True`,  the variants to annotate will be selected automatically using a sliding window with `windowsize=500`kb. chr:pos
    * `string`: the column name used for annotation
    * `"GENENAME"` : automatically annotate nrearest gene names, using pyensembl. (remember to specify `build`, default is `build="19"`)
* `repel_force` : when the annotation overlaps with other, try increasing the repel_force to increase the padding between annotations. (deault is 0.01)
* `anno_set `: if you want to annotate only a few specific variants, you can simply provide a list of SNPIDs. 
* `anno_alias` : snpid:text dictionary for customized annotation

<img width="376" alt="image" src="https://user-images.githubusercontent.com/40289485/196592902-36c52102-09f2-4c58-894b-a7177e87bd09.png">


#### example:
`mysumstats.plot_mqq(skip=3,anno=True)`


<img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/196591270-592ce820-c541-4b84-a766-0b58bc6423ee.png">


`mysumstats.plot_mqq(skip=3,anno="GENENAME",build="19")`

<img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/196591371-262d31d5-9640-474f-af0d-d6c511c77280.png">

`mysumstats.plot_mqq(skip=3, anno_set=["rs12509595","19:15040733:T:C"])`

<img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/196591966-c9618c45-456b-4eb8-991b-66420f847a97.png">

`mysumstats.plot_mqq(skip=3, anno_set=["rs12509595","19:15040733:T:C"], anno_alias={"rs12509595":"anything you want here"})`

<img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/196592136-f2cbc488-e02f-409b-b0ac-2ab16f7b1fd4.png">


----------------------------------------------------------------

### Highlight specified loci

Highlight specified loci. 
- `highlight ` : specify the lead variants of loci for highlighting.
- `highlight_color`: specify the color ussed for highlighting.
- `highlight_windowkb` : specify the span of highlighted region (deault: `highlight_windowkb = 500` kp)

### Pinpoint specified variants

Pinpoint certain variants in the manhattan plot.
- `pinpoint` : a list of SNPIDs
- `pinpoint_color` : color for pinpoint

```
mysumstats.plot_mqq(skip=3,anno="GENENAME",build="19",
                   highlight=["rs12509595","rs7989823"],
                   pinpoint=["rs671","19:15040733:T:C"])
```
<img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/196593330-1794223c-79fd-40f9-942d-8acbbe00827b.png">

----------------------------------------------------------------

### Maf-stratified QQ plot
- `stratified`: if True, plot MAF straitified QQ plot. Require EAF in sumstats.
- `maf_bins`: maf bins for straitification. (default: maf_bins=[(0, 0.01), (0.01, 0.05), (0.05, 0.25),(0.25,0.5)])
- `maf_bin_colors`: colors used for each bin. (default: maf_bin_colors = ["#f0ad4e","#5cb85c", "#5bc0de","#000042"])

<img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/196593456-93539e3e-e8c0-476f-ae4e-e5843f970ba6.png">

----------------------------------------------------------------

### Use rank or POS

`use_rank`: if True, use the rank instead of real base pair position.

use simply rank and basepair position to draw the x axis.

### Colors and fonts

----------------------------------------------------------------


