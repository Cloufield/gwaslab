# Batch load LDSC log file

Simply batch load LDSc data into a pandas DataFrame for other manipulation.

GWASLab uses regular expression to match the values and fill them into a dataframe.

`gl.read_ldsc()` gl.read_ldsc(filelist, mode="h2")

```
`filelist` : a list of paths to ldsc log files

`mode` : h2 or rg 

```python
#mode=h2
ldsc_file_list=["file1."]
myldscrg = gl.read_ldsc(ldsc_file_list, mode=h2)


#mode=rg
ldsc_file_list=["file1."]
myldscrg = gl.read_ldsc(mode=h2)


pathlist=["./test.results.log","./test2.results.log"]

ldsc_h2 = gl.read_ldsc(pathlist, mode="h2")
ldsc_rg = gl.read_ldsc(pathlist, mode="rg")

ldsc_h2
Filename	h2_obs	h2_se	Lambda_gc	Mean_chi2	Intercept	Intercept_se	Ratio	Ratio_se
test.results.log	42.9954	8.657	1.2899	1.3226	0.0098	0.0098	0.6538	0.0304
test2.results.log	NA	NA	1.2899	1.3226	0.0098	0.0098	Ratio < 0	NA

ldsc_rg
p1	p2	rg	se	z	p	h2_obs	h2_obs_se	h2_int	h2_int_se	gcov_int	gcov_int_se
./test.results.log	./test.results.log	0.2317	0.0897	2.5824	0.0098	0.3305	0.0571	0.9612	0.009	-0.0001	0.0062
./test.results.log	./test2.results.log	0.2317	0.0897	2.5824	0.0098	0.3305	0.0571	0.9612	0.009	-0.0001	0.0062

```

For genetic correlation, after loading, you can use `gl.plot_rg()` to plot a heat map to visualize the results. No extra manipulation needed.

```python
gl.plot_rg(myldscrg)
```
