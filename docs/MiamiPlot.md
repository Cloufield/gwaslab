# Miami plot

As a standalone function, gwaslab can plot miami plot given a pair of sumstats files.

See examples [here.](https://cloufield.github.io/gwaslab/visualization_miami/)

## Usage

```
gl.plot_miami( 
          path1,
          path2,
          cols1=["CHR","POS","P"],
          cols2=["CHR","POS","P"],
          sep=["\t","\t"],
          anno= None
          )
```

- `path1` : `string` path to sumstats1
- `path2` : `string` path to sumstats2
- `cols1` : `list` CHR,POS,P names for sumstats1
- `cols2` : `list` CHR,POS,P names for sumstats2
- `sep`   : `list` separator for each sumstats
- `anno`  : `boolean` or `GENENAME`
- `region` : only plot a region. For example, `region=(2,2153874,21753874)` .
- `highlight1`: `list` , list of loci to highlight for sumstats1. For example, `highlight1=[(5,124289158)]` .
- `highlight2`: `list` , list of loci to highlight for sumstats2. 
- `pinpoint1`: `list` , list of variants to pinpolint for sumstats1.
- `pinpoint2`: `list`,  list of variants to pinpolint for sumstats2.
- `titles`:`list`, two titles. For example,`titles=["male","female"]`

## Example

```
mysumstats = gl.plot_miami(path1="bmi_male_bbj.txt.gz" ,
                           path2="bmi_female_bbj.txt.gz",
                           cols1=["CHR","POS","P"],
                           cols2=["CHR","POS","P"],
                           titles=["bmi male","bmi female"],
                           titles_pad=[0.15,0.0],
                           anno="GENENAME",
                           region_grid=True,
                           highlight1=[(5,124289158)],
                           pinpoint2=[(2,653874)]
                           )
```

<img width="700" alt="image" src="https://user-images.githubusercontent.com/40289485/197526569-7850041d-e247-4f69-8505-ef7750a6d4de.png">

```
mysumstats = gl.plot_miami(path1="bmi_male_bbj.txt.gz" ,
                           path2="bmi_female_bbj.txt.gz",
                           cols1=["CHR","POS","P"],
                           cols2=["CHR","POS","P"],
                           titles=["male","female"],
                           region=(2,2153874,21753874),
                           titles_pad=[0.1,-0.05],
                           anno="GENENAME",
                           region_grid=True
                           )
```

<img width="700" alt="image" src="https://user-images.githubusercontent.com/40289485/197531327-8fdd51d0-e21f-4da2-9e2d-c4f8a8bd5285.png">

Note: implemented in v3.3.4
