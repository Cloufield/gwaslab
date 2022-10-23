# Miami plot (not finished)

```
gl.plot_miami( path1="./t2d_bbj.txt.gz" ,
              path2="./hum0197.v3.BBJ.IS.v1/GWASsummary_IS_Japanese_SakaueKanai2020.auto.txt.gz" ,
              cols1=["CHR","POS","P"],
              cols2=["CHR","POS","p.value"],
              sep=["\t","\s+"],
              skip=1,cut=True)
```

- `path1` : `string` path to sumstats1
- `path2` : `string` path to sumstats2
- `cols1` : `list` CHR,POS,P names for sumstats1
- `cols2` : `list` CHR,POS,P names for sumstats2
- `sep`   : `list` separator for each sumstats

<img width=700 src="https://user-images.githubusercontent.com/40289485/197393700-901c5cec-92a4-4aa7-b299-a760b76c3108.png">
