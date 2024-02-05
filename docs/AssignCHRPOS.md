# Assigning CHR/POS

GWASLab can update CHR/POS using a pre-processed SNPID-rsID table to assign CHR and POS based on rsID. Currently, only variants in 1KG phase3v5 are supported.  

## .rsid_to_chrpos()

### Reference data
GWASLab provides a download function `gl.download_ref()` and two curated tables which contains ~80M 1KG variants:

- `hg19` : `gl.download_ref("1kg_dbsnp151_hg19_auto")`
- `hg38` : `gl.download_ref("1kg_dbsnp151_hg38_auto")`

[Reference data](https://cloufield.github.io/gwaslab/Reference/)

!!! note "1kg_dbsnp151_hg19_auto"
    ```
    ~/.gwaslab$ zcat 1kg_dbsnp151_hg19_auto.txt.gz |head
    SNPID   rsID    CHR     POS     NEA     EA
    1:10177:A:AC    rs367896724     1       10177   A       AC
    1:10235:T:TA    rs540431307     1       10235   T       TA
    1:10352:T:TA    rs555500075     1       10352   T       TA
    1:10505:A:T     rs548419688     1       10505   A       T
    1:10511:G:A     rs534229142     1       10511   G       A
    1:10539:C:A     rs537182016     1       10539   C       A
    1:10542:C:T     rs572818783     1       10542   C       T
    1:10579:C:A     rs538322974     1       10579   C       A
    1:10616:CCGCCGTTGCAAAGGCGCGCCG:C        rs376342519     1       10616   CCGCCGTTGCAAAGGCGCGCCG  C
    ```


### .rsid_to_chrpos()

```

# if not downloaded yet :
# gl.download_ref("1kg_dbsnp151_hg19_auto")

# assign chr and pos using rsID
mysumstats.rsid_to_chrpos( path = gl.get_path("1kg_dbsnp151_hg19_auto"))
```

| Option | DataType | Description                        | Default |
|--------|----------|------------------------------------|---------|
| `path` | `string` | the path to reference tabular file | -       |


### Example 

!!! example "Assign CHR and POS using rsID"

    ```
    mysumstats.data
    
    EA	NEA	EAF	BETA	SE	P	N	DIRECTION	STATUS	rsID
    0	G	A	0.9960	-0.0737	0.1394	0.59700	166718	-?+-	9960099	rs565766235
    1	G	A	0.0040	0.0737	0.1394	0.59730	166718	+?-+	9960099	rs534711480
    2	C	T	0.0051	0.0490	0.1231	0.69080	166718	+?-+	9960099	rs540210562
    3	TAA	T	0.8374	0.0213	0.0199	0.28460	166718	-?++	9960399	rs529266287
    4	T	A	0.8593	0.0172	0.0156	0.27050	166718	-?++	9960099	rs28544273
    ...	...	...	...	...	...	...	...	...	...	...
    9995	C	T	0.1292	-0.0350	0.0191	0.06686	191764	----	9960099	rs55938238
    9996	C	T	0.4886	-0.0014	0.0094	0.88070	191764	-0+-	9960099	rs7520225
    9997	C	T	0.9476	-0.0061	0.0216	0.77790	191764	---+	9960099	rs151238770
    9998	C	CT	0.9418	-0.0047	0.0199	0.81500	191764	---+	9960399	rs137909285
    9999	G	A	0.9828	-0.0084	0.0433	0.84620	191764	+--+	9960099	rs77575110
    
    
    mysumstats.rsid_to_chrpos(path=gl.get_path("1kg_dbsnp151_hg19_auto"))
    
    Fri Jan 27 00:06:24 2023 Start to update chromosome and position information based on rsID...
    Fri Jan 27 00:06:24 2023  -Current Dataframe shape : 10000  x  10
    Fri Jan 27 00:06:24 2023  -rsID dictionary file: /home/he/.gwaslab/1kg_dbsnp151_hg19_auto.txt.gz
    Fri Jan 27 00:06:24 2023  -Setting block size:  10000000
    Fri Jan 27 00:06:24 2023  -Loading block: 0   1   2   3   4   5   6   7   
    Fri Jan 27 00:11:44 2023  -Updating CHR and POS finished.Start to re-fixing CHR and POS... 
    Fri Jan 27 00:11:44 2023 Start to fix chromosome notation...
    Fri Jan 27 00:11:44 2023  -Current Dataframe shape : 10000  x  12
    Fri Jan 27 00:11:45 2023  -Variants with standardized chromosome notation: 9942
    Fri Jan 27 00:11:45 2023  -Variants with fixable chromosome notations: 0
    Fri Jan 27 00:11:45 2023  -Variants with NA chromosome notations: 58
    Fri Jan 27 00:11:45 2023  -No unrecognized chromosome notations...
    Fri Jan 27 00:11:46 2023 Finished fixing chromosome notation successfully!
    Fri Jan 27 00:11:46 2023 Start to fix basepair positions...
    Fri Jan 27 00:11:46 2023  -Current Dataframe shape : 10000  x  12
    Fri Jan 27 00:11:46 2023  -Converting to Int64 data type ...
    Fri Jan 27 00:11:46 2023  -Position upper_bound is: 250,000,000
    Fri Jan 27 00:11:46 2023  -Remove outliers: 0
    Fri Jan 27 00:11:46 2023  -Converted all position to datatype Int64.
    Fri Jan 27 00:11:46 2023 Finished fixing basepair position successfully!
    
    mysumstats.data
    
    rsID	EA	NEA	EAF	BETA	SE	P	N	DIRECTION	STATUS	CHR	POS
    0	rs565766235	G	A	0.9960	-0.0737	0.1394	0.59700	166718	-?+-	9960099	1	725932
    1	rs534711480	G	A	0.0040	0.0737	0.1394	0.59730	166718	+?-+	9960099	1	725933
    2	rs540210562	C	T	0.0051	0.0490	0.1231	0.69080	166718	+?-+	9960099	1	737801
    3	rs529266287	TAA	T	0.8374	0.0213	0.0199	0.28460	166718	-?++	9960399	1	749963
    4	rs28544273	T	A	0.8593	0.0172	0.0156	0.27050	166718	-?++	9960099	1	751343
    ...	...	...	...	...	...	...	...	...	...	...	...	...
    9995	rs55938238	C	T	0.1292	-0.0350	0.0191	0.06686	191764	----	9960099	1	3142135
    9996	rs7520225	C	T	0.4886	-0.0014	0.0094	0.88070	191764	-0+-	9960099	1	3142137
    9997	rs151238770	C	T	0.9476	-0.0061	0.0216	0.77790	191764	---+	9960099	1	3142161
    9998	rs137909285	C	CT	0.9418	-0.0047	0.0199	0.81500	191764	---+	9960399	1	3142212
    9999	rs77575110	G	A	0.9828	-0.0084	0.0433	0.84620	191764	+--+	9960099	1	3142762
    ```


## .rsid_to_chrpos2()

**Available since v3.4.31**

`.rsid_to_chrpos2()` is a function to assign CHR and POS based on rsID using a HDF5 file derived from dbSNP reference VCF files.

### Refernce VCF from dbSNP

First, download reference VCF file from dbSNP ftp site.

```
# For example, dbSNP v155 hg19
GCF_000001405.25.gz (24G)
GCF_000001405.25.gz.tbi (2.9M)
```


### gl.process_ref_vcf()
```
gl.process_ref_vcf()
```

Process the vcf file and convert it to HDF5 file using `.process_ref_vcf()`. This step may take up to one or two hours.

| Option      | DataType | Description                                                 | Default              |
|-------------|----------|-------------------------------------------------------------|----------------------|
| `vcf`       | `string` | the path to dbSNP VCF file                                  | -                    |
| `directory` | `string` | the directory where you want output the converted HDF5 file | the same as VCF file |

```
directory="/home/yunye/work/gwaslab/examples/vcf_hd5/"
vcf = "/home/yunye/CommonData/Reference/ncbi_dbsnp/ncbi_dbsnp/db155/GCF_000001405.25.gz"

gl.process_ref_vcf(vcf=vcf,
                   directory=directory,
                   chr_dict=gl.get_NC_to_number(build="19"))

Fri Nov 10 11:27:59 2023 Start processing VCF files:
Fri Nov 10 11:27:59 2023  -Reference VCF path:/home/yunye/CommonData/Reference/ncbi_dbsnp/ncbi_dbsnp/db155/GCF_000001405.25.gz
Fri Nov 10 11:27:59 2023  -Output group size:20000000
Fri Nov 10 11:27:59 2023  -Compression level:9
Fri Nov 10 11:27:59 2023  -Loading chunksize:20000000
Fri Nov 10 11:27:59 2023  -HDF5 Output path: /home/yunye/work/gwaslab/examples/vcf_hd5/rsID_CHR_POS_groups_20000000.h5
Fri Nov 10 11:27:59 2023  -Log output path: /home/yunye/work/gwaslab/examples/vcf_hd5/rsID_CHR_POS_groups_20000000.log
Fri Nov 10 11:27:59 2023  -Processing chunk: 0 1 2 3 4 ...
```

Assign CHR POS using the HDF5 file. 

### .rsid_to_chrpos2()

```
.rsid_to_chrpos2()
```


### Options
| Option    | DataType | Description                          | Default |
|-----------|----------|--------------------------------------|---------|
| `path`    | `string` | the path to the HDF5 file            | -       |
| `n_cores` | `int`    | number of threads to use             | `./`    |
| `build`   | `string` | genome build version for CHR and POS | `"99"`  |


### Example

```
...
mysumstats.rsid_to_chrpos2(path="/home/yunye/work/gwaslab/examples/vcf_hd5/rsID_CHR_POS_groups_20000000.h5")

Fri Nov 10 17:30:44 2023 Start to assign CHR and POS using rsIDs... 
Fri Nov 10 17:30:44 2023  -Source hdf5 file:  ./vcf_hd5/rsID_CHR_POS_groups_20000000.hd5
Fri Nov 10 17:30:44 2023  -Cores to use :  4
Fri Nov 10 17:30:44 2023  -Blocksize (make sure it is the same as hdf5 file ):  20000000
Fri Nov 10 17:30:44 2023  -Non-Valid rsIDs:  58
Fri Nov 10 17:30:44 2023  -Duplicated rsIDs except for the first occurrence:  0
Fri Nov 10 17:30:44 2023  -Valid rsIDs:  9942
Fri Nov 10 17:30:44 2023  -Initiating CHR ... 
Fri Nov 10 17:30:44 2023  -Initiating POS ... 
Fri Nov 10 17:30:44 2023  -Divided into groups:  41
Fri Nov 10 17:30:44 2023   - {0, 1, 2, 3, 4, 5, 6, 7, 9, 10, 18, 19, 26, 27, 28, 37, 38, 39, 43, 44, 48, 49, 52, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74}
Fri Nov 10 17:30:46 2023  -Number of groups in HDF5:  92
Fri Nov 10 17:30:46 2023  -Max index of groups in HDF5:  105
Fri Nov 10 17:32:44 2023  -Merging group data... 
Fri Nov 10 17:32:44 2023  -Append data... 
Fri Nov 10 17:32:44 2023 Start to fix chromosome notation...
Fri Nov 10 17:32:44 2023  -Current Dataframe shape : 10000  x  13
Fri Nov 10 17:32:44 2023  -Checking CHR data type...
Fri Nov 10 17:32:44 2023  -Variants with standardized chromosome notation: 9705
Fri Nov 10 17:32:44 2023  -Variants with fixable chromosome notations: 0
Fri Nov 10 17:32:44 2023  -Variants with NA chromosome notations: 295
Fri Nov 10 17:32:44 2023  -No unrecognized chromosome notations...
Fri Nov 10 17:32:44 2023  -Sanity check for CHR...
Fri Nov 10 17:32:44 2023  -Removed 0 variants with CHR < 1...
Fri Nov 10 17:32:44 2023 Finished fixing chromosome notation successfully!
Fri Nov 10 17:32:44 2023 Start to fix basepair positions...
Fri Nov 10 17:32:44 2023  -Current Dataframe shape : 10000  x  13
Fri Nov 10 17:32:44 2023  -Converting to Int64 data type ...
Fri Nov 10 17:32:44 2023  -Position upper_bound is: 250,000,000
Fri Nov 10 17:32:44 2023  -Remove outliers: 0
Fri Nov 10 17:32:44 2023  -Converted all position to datatype Int64.
Fri Nov 10 17:32:44 2023 Finished fixing basepair position successfully!
Fri Nov 10 17:32:44 2023 Finished assigning CHR and POS using rsIDs.
```