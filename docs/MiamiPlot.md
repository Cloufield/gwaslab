# Miami plot (not finished)

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



# Example

```
mysumstats = gl.plot_miami(path1="bmi_male_bbj.txt.gz" ,
                           path2="bmi_female_bbj.txt.gz",
                           cols1=["CHR","POS","P"],
                           cols2=["CHR","POS","P"],
                           titles=["male","female"],
                           titles_pad=[0.15,0.0],
                           readcsv_args={"nrows":100000},
                           anno="GENENAME",
                           highlight=[(5,124289158)],
                           pinpoint=[(2,653874)]
                           )
```

<img width="700" alt="image" src="https://user-images.githubusercontent.com/40289485/197463243-89352749-f882-418d-907d-27530fd4e922.png">

