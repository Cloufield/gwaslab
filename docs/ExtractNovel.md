
# Checking if lead variants are novel or not

GWASLab can check if the lead variants of your summary statistics overlap with reported variants or not based on the physical distance.

## get_novel()

```
sumstats.getnovel(
                    known,
                    efo,
                    only_novel=False,
                    windowsizekb_for_novel=1000,
                    windowsizekb=500,
                    sig_level=5e-8,
                    output_known=False)
```

Required (either):

- `known` : `string`, path to local file of reported variants
- `efo` : `string`, efo id for the target trait, which is used for querying the GWASCatalog.

Options:

- `sig_level=5e-8` : significance threshold for extracting lead variants
- `windowsizekb=500` : windowsize for extracting lead variants
- `windowsizekb_for_novel=1000` : windowsize for determining if lead variants overlap with reported variants in GWASCatalog
- `only_novel`: output only novel variants
- `output_known`: additionally output the reported variants



## Checking with GWAScatalog API 

!!! warning "Only works when your sumstats are based on GRCh38" 

!!! example "Querying the GWAS Catalog"

    ```python
    # sample data
    Reference : Suzuki, K., Akiyama, M., Ishigaki, K., Kanai, M., Hosoe, J., Shojima, N., ... & Kadowaki, T. (2019). Identification of 28 new     susceptibility loci for type 2 diabetes in the Japanese population. Nature genetics, 51(3), 379-386.
    !wget -O t2d_bbj.txt.gz http://jenger.riken.jp/14/
    
    # load data and run a basic check
    mysumstats = gl.Sumstats("t2d_bbj.txt.gz",
                 snpid="SNP",
                 chrom="CHR",
                 pos="POS",
                 se="SE",
                 p="P")
    mysumstats.basic_check()

    # liftover to grch38 (original sumstats are based on grch37)
    mysumstats.liftover(from_build="19",to_build="38")

    # run get_novel with efo id for type 2 diabetes mellitus
    mysumstats.get_novel(efo="MONDO_0005148")

    Fri Feb  3 12:11:44 2023 Start to check if lead variants are known...
    Fri Feb  3 12:11:44 2023 Start to extract lead variants...
    Fri Feb  3 12:11:44 2023  -Processing 1035234 variants...
    Fri Feb  3 12:11:44 2023  -Significance threshold : 5e-08
    Fri Feb  3 12:11:44 2023  -Sliding window size: 500  kb
    Fri Feb  3 12:11:44 2023  -Found 9458 significant variants in total...
    Fri Feb  3 12:11:45 2023  -Identified 89 lead variants!
    Fri Feb  3 12:11:45 2023 Finished extracting lead variants successfully!
    Fri Feb  3 12:11:45 2023 Start to retrieve data from GWASCatalog...
    Fri Feb  3 12:11:45 2023  -Requesting (GET) trait information through the GWASCatalog API...
    Fri Feb  3 12:11:45 2023  -EFO trait api: https://www.ebi.ac.uk/gwas/rest/api/efoTraits/MONDO_0005148
    Fri Feb  3 12:11:46 2023  -Trait Name: type 2 diabetes mellitus
    Fri Feb  3 12:11:46 2023  -Trait URL: http://purl.obolibrary.org/obo/MONDO_0005148
    Fri Feb  3 12:11:46 2023  -Requesting (GET) GWAS associations through the GWASCatalog API...
    Fri Feb  3 12:11:46 2023  -associationsByTraitSummary API: https://www.ebi.ac.uk/gwas/rest/api/efoTraits/MONDO_0005148/associations?    projection=associationByEfoTrait
    Fri Feb  3 12:11:46 2023  -Note: this step might take a while...
    Fri Feb  3 12:18:55 2023  -Status code 200 OK: Retrieved data from GWASCatalog successffully ...
    Fri Feb  3 12:18:55 2023  -Loading json ...
    Fri Feb  3 12:18:58 2023  -Parsing json ...
    Fri Feb  3 12:18:58 2023  -Number of reported associations for MONDO_0005148 in GWASCatalog: 5345
    Fri Feb  3 12:18:58 2023  -Loading retrieved data into gwaslab Sumstats object ...
    Fri Feb  3 12:18:58 2023 GWASLab version 3.3.24 https://cloufield.github.io/gwaslab/
    Fri Feb  3 12:18:58 2023 (C) 2022-2023, Yunye He, Kamatani Lab, MIT License, gwaslab@gmail.com
    Fri Feb  3 12:18:58 2023   - format_name  : gwaslab
    Fri Feb  3 12:18:58 2023   - format_source  : https://cloufield.github.io/gwaslab/
    Fri Feb  3 12:18:58 2023   - format_version  : 20220729_v3
    Fri Feb  3 12:19:00 2023 Finished retrieving data from GWASCatalog...
    Fri Feb  3 12:19:00 2023  -Retrieved 4036 associations from GWAS catalog.
    Fri Feb  3 12:19:00 2023  -Lead variants in known loci: 4036
    Fri Feb  3 12:19:00 2023  -Checking the minimum distance between identified lead variants and provided known variants...
    Fri Feb  3 12:19:00 2023  -Identified  89  known vairants in current sumstats...
    Fri Feb  3 12:19:00 2023  -Identified  0  novel vairants in current sumstats...
    Fri Feb  3 12:19:00 2023 Finished checking known or novel successfully!
    
    SNPID	CHR	POS	SE	P	MLOG10P	STATUS	DISTANCE_TO_KNOWN	KNOWN_ID	KNOWN_PUBMED_ID	KNOWN_AUTHOR	NOVEL	LOCATION_OF_KNOWN
    0	1:22068326_A_G	1	21741833	0.0103	1.629000e-09	8.788079	3860999	0	rs1825307	30718926	Suzuki K	False	Same
    1	1:51103268_T_C	1	50637596	0.0120	2.519000e-11	10.598772	3860999	0	rs12031188	30718926	Suzuki K	False	Same
    2	1:154309595_TA_T	1	154337119	0.0166	3.289000e-08	7.482936	3860999	1	rs68062313	30718926	Suzuki K	False	    Upstream
    3	2:640986_CACAT_C	2	640986	0.0150	2.665000e-10	9.574303	3860999	1	rs72156956	30718926	Suzuki K	False	Upstream
    4	2:27734972_G_A	2	27512105	0.0088	3.897000e-15	14.409270	3860999	0	rs6547692	30718926	Suzuki K	False	Same
    ...	...	...	...	...	...	...	...	...	...	...	...	...	...
    84	X:21569920_A_G	23	21551802	0.0076	2.616000e-08	7.582362	3860999	0	rs6633421	30718926	Suzuki K	False	Same
    85	X:48724648_CAA_C	23	48866248	0.0103	4.576000e-09	8.339514	3860999	1	rs782100977	30718926	Suzuki K	False	    Upstream
    86	X:57170781_A_AT	23	57144348	0.0076	4.583000e-09	8.338850	3860999	1	rs144226500	30718926	Suzuki K	False	Upstream
    87	X:117915163_T_TA	23	118781200	0.0071	9.818000e-15	14.007977	3860999	1	rs11390176	30718926	Suzuki K	False	    Upstream
    88	X:152908887_G_A	23	153643433	0.0077	9.197000e-58	57.036354	3860999	0	rs1894299	30718926	Suzuki K	False	Same
    ```
    
    All of the lead variants are reported which is expected of course.

!!! warning 
    GWAS Catalog API is unstable sometimes.

## Checking with local files

!!! example
    
    Suppose we have a list of reported variants like:

    ```bash
    cat ./toy_data/known_loci.txt
    CHR POS
    1 154309595
    1 51103268
    ```

    ```python
    mysumstats.get_novel(known="./toy_data/known_loci.txt")
    
    Fri Feb  3 11:44:12 2023 Start to check if lead variants are known...
    Fri Feb  3 11:44:12 2023 Start to extract lead variants...
    Fri Feb  3 11:44:12 2023  -Processing 1035804 variants...
    Fri Feb  3 11:44:12 2023  -Significance threshold : 5e-08
    Fri Feb  3 11:44:12 2023  -Sliding window size: 500  kb
    Fri Feb  3 11:44:12 2023  -Found 9461 significant variants in total...
    Fri Feb  3 11:44:13 2023  -Identified 89 lead variants!
    Fri Feb  3 11:44:13 2023 Finished extracting lead variants successfully!
    Fri Feb  3 11:44:13 2023  -Lead variants in known loci: 2
    Fri Feb  3 11:44:13 2023  -Checking the minimum distance between identified lead variants and provided known variants...
    Fri Feb  3 11:44:13 2023  -Identified  2  known vairants in current sumstats...
    Fri Feb  3 11:44:13 2023  -Identified  87  novel vairants in current sumstats...
    Fri Feb  3 11:44:13 2023 Finished checking known or novel successfully!
    
    SNPID	CHR	POS	SE	P	MLOG10P	STATUS	DISTANCE_TO_KNOWN	KNOWN_ID	NOVEL	LOCATION_OF_KNOWN
    0	1:22068326_A_G	1	22068326	0.0103	1.629000e-09	8.788079	9960999	29034942	1:51103268	True	Upstream
    1	1:51103268_T_C	1	51103268	0.0120	2.519000e-11	10.598772	9960999	0	1:51103268	False	Same
    2	1:154309595_TA_T	1	154309595	0.0166	3.289000e-08	7.482936	9960999	0	1:154309595	False	Same
    3	2:640986_CACAT_C	2	640986	0.0150	2.665000e-10	9.574303	9960999	<NA>	<NA>	True	NoneOnThisChr
    4	2:27734972_G_A	2	27734972	0.0088	3.897000e-15	14.409270	9960999	<NA>	<NA>	True	NoneOnThisChr
    ...	...	...	...	...	...	...	...	...	...	...	...
    84	X:21569920_A_G	23	21569920	0.0076	2.616000e-08	7.582362	9960999	<NA>	<NA>	True	NoneOnThisChr
    85	X:48724648_CAA_C	23	48724648	0.0103	4.576000e-09	8.339514	9960999	<NA>	<NA>	True	NoneOnThisChr
    86	X:57170781_A_AT	23	57170781	0.0076	4.583000e-09	8.338850	9960999	<NA>	<NA>	True	NoneOnThisChr
    87	X:117915163_T_TA	23	117915163	0.0071	9.818000e-15	14.007977	9960999	<NA>	<NA>	True	NoneOnThisChr
    88	X:152908887_G_A	23	152908887	0.0077	9.197000e-58	57.036354	9960999	<NA>	<NA>	True	NoneOnThisChr
    ```
