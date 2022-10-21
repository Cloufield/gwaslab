# Regional plot:

<img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/197126045-b1c55adf-3391-4c3d-b2f6-eaeac7c26024.png">

```
mqqplot(  vcf_path=None,
          vcf_chr_dict=get_number_to_chr(),
          gtf_path="defualt",
          mlog10p="MLOG10P",
          scaled=False,
          mode="r",
          # region
          region = None,
          region_step = 21,
          region_grid = False,
          region_grid_line = {"linewidth": 2,"linestyle":"--"},
          region_lead_grid = True,
          region_lead_grid_line = {"alpha":0.5,"linewidth" : 2,"linestyle":"--","color":"#FF0000"},
          region_hspace=0.02,
          region_ld_threshold = [0.2,0.4,0.6,0.8],
          region_ld_colors = ["#E4E4E4","#020080","#86CEF9","#24FF02","#FDA400","#FF0000","#FF0000"],
          region_recombination = True,
          region_protein_coding=True,
          region_flank_factor = 0.05,
          taf=[4,0,0.95,1,1],
          # track_n, track_n_offset,font_ratio,exon_ratio,text_offset
          mqqratio=3,
          windowsizekb=500,
          anno=None,
          anno_set=[],
          anno_alias={},
          anno_d={},
          anno_source = "ensembl",
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
          log=Log()
          ):
```
gwaslab regional plot function is based on plot_mqq().

## Regional mqq plot
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


## Regional plot without LD
```
mysumstats.plot_mqq(mode="r", region=(7,156538803,157538803),region_grid=True)
```
<img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/197126849-72439c54-19fe-4fb9-a795-e35f42a9236b.png">



## Regional plot
```
mysumstats.plot_mqq(mode="r",region=(7,156538803,157538803),region_grid=True,
                    vcf_path="/home/yunye/mydata/d_disk/eas_1kg_af/EAS.chr7.split_norm_af.vcf.gz")
```
<img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/197126045-b1c55adf-3391-4c3d-b2f6-eaeac7c26024.png">
