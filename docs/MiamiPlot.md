# Miami plot

!!! info "Implemented since v3.3.4"

As a standalone function, GWASLab can plot miami plot given a pair of sumstats files.

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

!!! example "Miami plot for male- amd female-specific GWAS on BMI"
    
    Sample datasets are obtained from JENGER:
    
    ```
    !wget -O bmi_male_bbj.txt.gz http://jenger.riken.jp/2analysisresult_qtl_download/
    !wget -O bmi_female_bbj.txt.gz http://jenger.riken.jp/4analysisresult_qtl_download/
    ```
    
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
