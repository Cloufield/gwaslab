# Reserved headers

Here is a list of all headers (gl.Sumstats.data) used in GWASLab functions. These headers are fixed to avoid ambiguity within GWASLab framework.

| Header | Keyword of args in preformat | Datatype | Description | Note |
|--------|----------------------------|----------|-------------|------|
| `SNPID` | snpid | `string` | variant ID (CHR:POS:NEA:EA) | GWASLab assumes EA=ALT and NEA=REF |
| `rsID` | rsid | `string` | dbSNP rsID | - |
| `CHR` | chrom | `Int64` | chromosome number (X 23, Y 24, MT 25) | - |
| `POS` | pos | `Int64` | base pair position | - |
| `EA` | ea | `category` | effect allele | - |
| `NEA` | nea | `category` | non-effect allele | - |
| `STATUS` | status | `category` | 7-digit variant status | lookup_status() to check |
| `REF` | ref | `category` | reference allele in reference genome | - |
| `ALT` | alt | `category` | alternative allele | - |
| `EAF` | eaf | `float64` | effect allele frequency | - |
| `RAF` | raf | `float64` | reference allele frequency (from reference VCF/BCF file) | typically annotated by infer_strand2() or annotate_sumstats() |
| `NEAF` | neaf | `float64` | non-effect allele frequency | - |
| `MAF` | maf | `float64` | minor allele frequency | if EAF, use fill_data to get MAF |
| `DAF` | daf | `float64` | difference in allele frequency (EAF - ALT_AF from reference VCF) | this DAF is not the derived allele frequency |
| `INFO` | info | `float32` | imputation INFO/RSQ | - |
| `BETA` | beta | `float64` | effect size beta | - |
| `SE` | se | `float64` | standard error of beta | - |
| `BETA_95U` | beta_95U | `float64` | upper bound of beta 95% condidence interval | - |
| `BETA_95L` | beta_95L | `float64` | lower bound of beta 95% condidence interval | - |
| `OR` | OR | `float64` | odds ratio | - |
| `OR_95U` | OR_95U | `float64` | upper bound of OR 95% condidence interval | - |
| `OR_95L` | OR_95L | `float64` | lower bound of OR 95% condidence interval | - |
| `HR` | HR | `float64` | hazard ratio | - |
| `HR_95U` | HR_95U | `float64` | upper bound of HR 95% condidence interval | - |
| `HR_95L` | HR_95L | `float64` | lower bound of HR 95% condidence interval | - |
| `CHISQ` | chisq | `float64` | chi square | - |
| `Z` | z | `float64` | z score | - |
| `T` | t | `float64` | t statistics | - |
| `F` | f | `float64` | F statistics | - |
| `P` | p | `float64` | P value | - |
| `P_MANTISSA` | p_mantissa | `float64` | P mantissa | - |
| `P_EXPONENT` | p_exponent | `float64` | P exponent | - |
| `MLOG10P` | mlog10p | `float64` | -log10(P) | - |
| `SNPR2` | snpr2 | `float64` | per variant R2 | - |
| `DOF` | dof | `Int64` | degree of freedom | - |
| `P_HET` | phet | `float64` | heterogeneity test P value | - |
| `I2_HET` | i2 | `float64` | heterogeneity I2 | - |
| `I2` | i2 | `float64` | heterogeneity I2 (alias for I2_HET) | - |
| `DIRECTION` | direction | `string` | direction of effect (+, -, ?) | - |
| `DENSITY` | density | `Int64` | signal density | - |
| `N` | n | `Int64` | total sample size | - |
| `N_CASE` | ncase | `Int64` | number of cases | - |
| `N_CONTROL` | ncontrol | `Int64` | number of controls | - |
| `N_EFF` | neff | `float64` | Effective sample size | - |
| `GENENAME` |   | `string` | nearest gene symbol | - |
| `CIS/TRANS` |   | `string` | whether the variant is in cis or trans region | Cis,Trans,NoReference |
| `DISTANCE_TO_KNOWN` |   | `Int64` | distance to nearest known variants | - |
| `LOCATION_OF_KNOWN` |   | `string` | relative location to nearest known variants | Same,Upstream,Downstream,NoReference |
| `KNOWN_ID` |   | `string` | nearest known variant ID | - |
| `KNOWN_PUBMED_ID` |   | `string` | pubmed ID of the known variant | - |
| `KNOWN_AUTHOR` |   | `string` | author of the study | - |
| `KNOWN_SET_VARIANT` |   | `string` | known set and overlapping variant | - |
| `KNOWN_VARIANT` |   | `string` | known variant overlapping with the variant | - |
| `KNOWN_SET` |   | `string` | variant set of the known variant | - |
| `NOVEL` |   | `string` | if the identified variants are novel | - |
| `PIP` | pip | `float64` | Posterior Inclusion Probability | - |
| `CREDIBLE_SET_INDEX` | credible_set_index | `Int64` | credible sets index | - |
| `N_SNP` | n_snp | `Int64` | number of variants included in this locus for finemapping | - |
| `LOCUS` | locus | `string` | locus name, usually the lead variant of the locus | - |
| `STUDY` | study | `string` | study name | - |
| `TEST` | test | `string` | test type | - |
| `BETA_RANDOM` | beta_random | `float64` | random effects beta | - |
| `SE_RANDOM` | se_random | `float64` | random effects standard error | - |
| `Z_RANDOM` | z_random | `float64` | random effects z score | - |
| `P_RANDOM` | p_random | `float64` | random effects P value | - |
