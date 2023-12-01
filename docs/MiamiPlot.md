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

|Option|DataType|Description|Default|
|-|-|-|-|
|`path1`|`string`, `pd.DataFrame`, or `gl.Sumstats`|path to sumstats1, or pd.DataFrame object, or gl.Sumstats Object|-|
|`path2`|`string`, `pd.DataFrame`, or `gl.Sumstats`|path to sumstats2, or pd.DataFrame object, or gl.Sumstats Object|-|
|`cols1`| `list`|CHR,POS,P names for sumstats1|`["CHR","POS","P"]`|
|`cols2`| `list`|CHR,POS,P names for sumstats2|`["CHR","POS","P"]`|
|`sep`|`list`|separator for each sumstats (when path is `string`)|`["\t","\t"]`|
|`anno`|`boolean` or `GENENAME`|if True, annotate `CHR:POS`; if `GENENAME`, annotate gene name|`None`|
|`region` |`tuple`|only plot a region. For example, `region=(2,2153874,21753874)` |`None`|
|`highlight`|`list`|list of loci (tuples of CHR and POS) to highlight for both sumstats. For example, `highlight=[(2,2153874),(5,124289158)]` |`None`|
|`highlight1`|`list`|list of loci (tuples of CHR and POS) to highlight for sumstats1. For example, `highlight1=[(2,2153874),(5,124289158)]`|`None`|
|`highlight2`|`list`|list of loci (tuples of CHR and POS) to highlight for sumstats2. For example, `highlight2=[(2,2153874),(5,124289158)]`|`None`|
|`highlight_windowkb`|`int`|window size for highlighting loci|`500`|
|`highlight_color`|`color`|color for highlighting loci|`"#CB132D"`|
|`pinpoint`|`list`|list of variants (tuples of CHR and POS) to pinpoint for both sumstats.  For example, `pinpoint=[(2,2153874),(5,124289158)]`|`None`|
|`pinpoint1`|`list`|list of variants (tuples of CHR and POS) to pinpoint for sumstats1.  For example, `pinpoint1=[(2,2153874),(5,124289158)]`|`None`|
|`pinpoint2`|`list`|list of variants (tuples of CHR and POS) to pinpoint for sumstats2.  For example, `pinpoint2=[(2,2153874),(5,124289158)]`|`None`|
|`pinpoint_color`|`color`|color for highlighting loci|`"red"`|
|`anno_set`|`list`|list of variants (tuples of CHR and POS) to annotate for both sumstats.  For example, `anno_set=[(2,2153874),(5,124289158)]`|`None`|
|`anno_set1`|`list`|list of variants (tuples of CHR and POS) to annotate for sumstats1.  For example, `anno_set1=[(2,2153874),(5,124289158)]`|`None`|
|`anno_set2`|`list`|list of variants (tuples of CHR and POS) to annotate for sumstats2.  For example, `anno_set2=[(2,2153874),(5,124289158)]`|`None`|
|`titles`|`list`|titles for sumstats. For example,`titles=["male","female"]`|`["",""]`|
|`titles_pad`|`list`|paddings for titles|`[0.2,0.2]`|
|`build`|`list`|genome build for annotate gene names. For example,`titles=["male","female"]`|`"19"`|
|`save`|`boolean` or `string`|if true, save to default path. if string, save to path specified by string|`True`|
|`save_args`|`dict`|additional parameters for `plt.save_fig`|`{"dpi":100,"facecolor":"white"}`|

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
                               highlight1=[(5,124289158)],
                               pinpoint2=[(2,653874)]
                               )
    ```
    
    <img width="700" alt="image" src="https://user-images.githubusercontent.com/40289485/197526569-7850041d-e247-4f69-8505-ef7750a6d4de.png">
    