# Brisbane plot

```python
import gwaslab as gl
import pandas as pd
```

```python
gl.show_version()
```

**stdout:**
```python
2024/12/23 12:19:09 GWASLab v3.5.4 https://cloufield.github.io/gwaslab/
2024/12/23 12:19:09 (C) 2022-2024, Yunye He, Kamatani Lab, MIT License, gwaslab@gmail.com
```

## Download data

Yengo, L., Vedantam, S., Marouli, E., Sidorenko, J., Bartell, E., Sakaue, S., ... & Lee, J. Y. (2022). A saturated map of common genetic variants associated with human height. Nature, 1-16.
Fig. 2: Brisbane plot

```python
!wget https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-022-05275-y/MediaObjects/41586_2022_5275_MOESM3_ESM.xlsx
```

**stdout:**
```python
--2024-12-23 12:19:13--  https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-022-05275-y/MediaObjects/41586_2022_5275_MOESM3_ESM.xlsx
Resolving static-content.springer.com (static-content.springer.com)... 146.75.112.95
Connecting to static-content.springer.com (static-content.springer.com)|146.75.112.95|:443... connected.
HTTP request sent, awaiting response... 200 OK
Length: 9368771 (8.9M) [application/octet-stream]
Saving to: ‘41586_2022_5275_MOESM3_ESM.xlsx.1’

41586_2022_5275_MOE 100%[===================>]   8.93M  46.0MB/s    in 0.2s    

2024-12-23 12:19:14 (46.0 MB/s) - ‘41586_2022_5275_MOESM3_ESM.xlsx.1’ saved [9368771/9368771]
```

## Read into pandas dataframe

```python
data = pd.read_excel("41586_2022_5275_MOESM3_ESM.xlsx",sheet_name=10,skiprows=1)
data
```

**stderr:**
```python
/home/yunye/anaconda3/envs/gwaslab_py39/lib/python3.9/site-packages/openpyxl/worksheet/_reader.py:329: UserWarning: Unknown extension is not supported and will be removed
  warn(msg)
```

```python
| Locus | Chr | SNP | BP (HG19) | BP (HG38) Effect Allele (A1) | \ |
| --- | --- | --- | --- | --- | --- |
| 0 | 1:METAFE | 1 | rs2710888 | 959842 | 1024462 |
| 1 | 1:METAFE | 1 | rs3934834 | 1005806 | 1070426 |
| 2 | 2:METAFE | 1 | rs182532 | 1287040 | 1351660 |
| 3 | 2:METAFE | 1 | rs17160669 | 1305561 | 1370181 |
| 4 | 3:METAFE | 1 | rs9660106 | 1797947 | 1866508 |
| ... | ... | ... | ... | ... | ... |
| 12106 | 0.764848 |  |  |  |  |
| 12107 | 0.831684 |  |  |  |  |
| 12108 | 0.379943 |  |  |  |  |
| 12109 | 0.751260 |  |  |  |  |
| 12110 | 0.935627 |  |  |  |  |

*[12111 rows x 26 columns]*
```

## Load into gwaslab Sumstats

```python
sumstats = gl.Sumstats(data,snpid="SNP",chrom="Chr",pos="BP (HG19)",p="b P-value")
```

**stdout:**
```python
2024/12/23 12:19:16 GWASLab v3.5.4 https://cloufield.github.io/gwaslab/
2024/12/23 12:19:16 (C) 2022-2024, Yunye He, Kamatani Lab, MIT License, gwaslab@gmail.com
2024/12/23 12:19:16 Start to initialize gl.Sumstats from pandas DataFrame ...
2024/12/23 12:19:16  -Reading columns          : BP (HG19),Chr,SNP,b P-value
2024/12/23 12:19:16  -Renaming columns to      : POS,CHR,SNPID,P
2024/12/23 12:19:16  -Current Dataframe shape : 12111  x  4
2024/12/23 12:19:16  -Initiating a status column: STATUS ...
2024/12/23 12:19:16  #WARNING! Version of genomic coordinates is unknown...
2024/12/23 12:19:17 Start to reorder the columns...v3.5.4
2024/12/23 12:19:17  -Current Dataframe shape : 12111 x 5 ; Memory usage: 21.88 MB
2024/12/23 12:19:17  -Reordering columns to    : SNPID,CHR,POS,P,STATUS
2024/12/23 12:19:17 Finished reordering the columns.
2024/12/23 12:19:17  -Column  : SNPID  CHR    POS   P       STATUS  
2024/12/23 12:19:17  -DType   : object string int64 float64 category
2024/12/23 12:19:17  -Verified: T      F      T     T       T       
2024/12/23 12:19:17  #WARNING! Columns with possibly incompatible dtypes: CHR
2024/12/23 12:19:17  -Current Dataframe memory usage: 21.88 MB
2024/12/23 12:19:17 Finished loading data successfully!
```

## Brisbane plot

- `mode=b` : Brisbane plot

```python
sumstats.plot_mqq(mode="b")
```

**stdout:**
```python
2024/12/23 12:19:17 Start to create MQQ plot...v3.5.4:
2024/12/23 12:19:17  -Genomic coordinates version: 99...
2024/12/23 12:19:17  #WARNING! Genomic coordinates version is unknown.
2024/12/23 12:19:17  -Genome-wide significance level to plot is set to 5e-08 ...
2024/12/23 12:19:17  -Raw input contains 12111 variants...
2024/12/23 12:19:17  -MQQ plot layout mode is : b
2024/12/23 12:19:17 Finished loading specified columns from the sumstats.
2024/12/23 12:19:17 Start data conversion and sanity check:
2024/12/23 12:19:17  -Removed 0 variants with nan in CHR or POS column ...
2024/12/23 12:19:17  -Removed 0 variants with CHR <=0...
2024/12/23 12:19:17  -Calculating DENSITY with windowsize of  100  kb
2024/12/23 12:19:17  -Converting data above cut line...
2024/12/23 12:19:17  -Maximum DENSITY value is 24.0 .
2024/12/23 12:19:17 Finished data conversion and sanity check.
2024/12/23 12:19:17 Start to create MQQ plot with 12111 variants...
2024/12/23 12:19:17  -Creating background plot...
2024/12/23 12:19:17 Finished creating MQQ plot successfully!
2024/12/23 12:19:17 Start to extract variants for annotation...
2024/12/23 12:19:17 Finished extracting variants for annotation...
2024/12/23 12:19:17 Start to process figure arts.
2024/12/23 12:19:17  -Processing X ticks...
2024/12/23 12:19:17  -Processing X labels...
2024/12/23 12:19:17  -Processing Y labels...
2024/12/23 12:19:17  -Processing Y tick lables...
2024/12/23 12:19:17  -Processing Y labels...
2024/12/23 12:19:17  -Processing lines...
2024/12/23 12:19:17  -Plotting horizontal line (  mean DENISTY): y = 1.9998348608702832
2024/12/23 12:19:17  -Plotting horizontal line ( median DENISTY): y = 1.0
2024/12/23 12:19:17 Finished processing figure arts.
2024/12/23 12:19:17 Start to annotate variants...
2024/12/23 12:19:17  -Skip annotating
2024/12/23 12:19:17 Finished annotating variants.
2024/12/23 12:19:17 Start to save figure...
2024/12/23 12:19:17  -Skip saving figure!
2024/12/23 12:19:17 Finished saving figure...
2024/12/23 12:19:17 Finished creating plot successfully!
```

```python
(<Figure size 3000x1000 with 1 Axes>, <gwaslab.g_Log.Log at 0x7f1446bc8b80>)
```

![Output image](images/notebooks/visualization_brisbane_img_0.png)

## Some customizations

```python
sumstats.plot_mqq(mode="b",anno=True,sig_line_color="red", windowsizekb=100000, verbose=False)
```

**stdout:**
```python
2024/12/23 12:19:18  #WARNING! Genomic coordinates version is unknown.
```

```python
(<Figure size 3000x1000 with 1 Axes>, <gwaslab.g_Log.Log at 0x7f1446bc8b80>)
```

![Output image](images/notebooks/visualization_brisbane_img_1.png)

## Annotate with gene names

```python
sumstats.plot_mqq(mode="b",anno="GENENAME",sig_line_color="red",build="19",windowsizekb=100000, verbose=False)
```

```python
(<Figure size 3000x1000 with 1 Axes>, <gwaslab.g_Log.Log at 0x7f1446bc8b80>)
```

![Output image](images/notebooks/visualization_brisbane_img_2.png)
