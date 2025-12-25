# Genetic correlation loading and plotting
## Heatmap: Genetic correlation matrix

Available since v3.4.15

```python
gl.plot_rg(ldsc)
```

!!! example "Simple plot"
    ```python
    # load ldsc log files
    ldsc_log_file_list = ["ldscrg.log"] 
    
    ldsc = gl.read_ldsc(ldsc_log_file_list, mode="rg")
    
    gl.plot_rg( ldsc )
    ```
    
    ![image](https://github.com/Cloufield/gwaslab/assets/40289485/46ea20c5-5e4d-4ad4-86ea-965f3a953518)

**Options:**

- `ldscrg` : `DataFrame`, results from ldsc-rg. 4 columns are required, `p1`,`p2`,`p`,`rg`.
- `p1`: `string`, column name for trait1, defaul: `p1`
- `p2`: `string`, column name for trait2 defaul: `p2`
- `rg`: `string`, column name for rg defaul: `rg`
- `p`: `string`, column name for p defaul: `p`
- `sig_levels`: `list`, default: `[0.05]`
- `panno`: `boolean`, default: `True`
- `corrections`: `list`,`non` no correction, `fdr` FDR, `bon` bonferroni , default: `["non","fdr","bom"]`
- `panno_texts`: `list`, text to annotate significant correlations, match the number of `corrections` times the number of `sig_levels`, default: `["*","**","***"]`
- `sort_key`: `function`, sort the columns , default: `None`
- `equal_aspect`: ``, defaul: `True`
- `fontsize`: ``, defaul: `10`
- `save`: `string` or `boolean`, defaul: `None`
- `save_args`: `dict`, defaul: `None`
- `full_cell`: `tuple`, threshold for full cell, default: `("fdr",0.05)`
- `yticklabel_args`: `dict`, default: `{"fontsize":10}`
- `xticklabel_args`: `dict`, default: `{"rotation":45,"horizontalalignment":"left", "verticalalignment":"bottom","fontsize":10}`
- `colorbar_args` :  `dict` , default:`{"shrink":0.82}`
- `cmap`: `cmap`, default:`matplotlib.cm.get_cmap('RdBu')`


!!! example "Customized plot"
    ```python
    # or load from tabular files
    ldsc = pd.read_csv("toy_data/input_rg.txt",sep="\t")
    ldsc
    #prepare the final dataset
    trait =  pd.read_csv("toy_data/trait_list.txt",sep="\t")
    trait["order"] = range(len(trait))
    order = trait["TRAIT"].values
    trait_set1 = trait.loc[trait["order"]>=59,"TRAIT"].values
    trait_set2 = trait.loc[trait["order"]<59,"TRAIT"].values
    ldsc = ldsc.loc[((ldsc["p1"].isin(trait_set1))&(ldsc["p2"].isin(trait_set2))) | ((ldsc["p1"].isin(trait_set2))&(ldsc["p2"].isin(trait_set1))),:]
    map_dic={order[i]:i+1 for i in range(len(order))}
    key=lambda x:x.map(map_dic)
    
    # plot
    df = gl.plot_rg( ldsc,
                sig_levels=[0.05],
                corrections =["non"],
                p="q",
                p1="p2",
                p2="p1",
                full_cell=("non",0.05),
                panno_texts=["*"],
                fig_args={"figsize":(15,15),"dpi":300},
                colorbar_args={"shrink":0.4},
                panno_args={"size":12,"c":"black"},fdr_method="i",
                fontsize=8,
                sort_key=key
                )
    ```
    
    ![image](https://github.com/Cloufield/gwaslab/assets/40289485/4f6c7dd8-5047-435a-8b08-5346bc168052)

    sample data source: https://github.com/mkanai/ldsc-corrplot-rg , Kanai, M., Akiyama, M., Takahashi, A., Matoba, N., Momozawa, Y., Ikeda, M., ... & Kamatani, Y. (2018). Genetic analysis of quantitative traits in the     Japanese population links cell types to complex human diseases. Nature genetics, 50(3), 390-400.
