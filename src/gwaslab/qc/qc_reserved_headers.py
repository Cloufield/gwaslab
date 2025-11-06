researved_header="""

Here is a list of all headers (gl.Sumstats.data) used in GWASLab functions. These headers are fixed to avoid ambiguity within GWASLab framework.

| Header              | Datatype   | Description                                    | Note                                         |
|---------------------|------------|------------------------------------------------|----------------------------------------------|
| `SNPID`             | `string`   | variant ID (CHR:POS:NEA:EA)                    | GWASLab assumes EA=ALT and NEA=REF           |
| `rsID`              | `string`   | dbSNP rsID                                     | -                                            |
| `CHR`               | `Int64`    | chromosome number (X 23, Y 24, MT 25)          | -                                            |
| `POS`               | `Int64`    | base pair position                             | -                                            |
| `EA`                | `category` | effect allele                                  | -                                            |
| `NEA`               | `category` | non-effect allele                              | -                                            |
| `STATUS`            | `category` | variant standardization & harmonization status | -                                            |
| `REF`               | `category` | reference allele in reference genome           | -                                            |
| `ALT`               | `category` | alternative allele                             | -                                            |
| `EAF`               | `float64`  | effect allele frequency                        | -                                            |
| `NEAF`              | `float64`  | non-effect allele frequency                    | -                                            |
| `MAF`               | `float64`  | minor allele frequency                         | -                                            |
| `INFO`              | `float32`  | imputation INFO/RSQ                            | -                                            |
| `BETA`              | `float64`  | effect size beta                               | -                                            |
| `SE`                | `float64`  | standard error of beta                         | -                                            |
| `BETA_95U`          | `float64`  | upper bound of beta 95% condidence interval    | -                                            |
| `BETA_95L`          | `float64`  | lower bound of beta 95% condidence interval    | -                                            |
| `OR`                | `float64`  | odds ratio                                     | -                                            |
| `OR_95U`            | `float64`  | upper bound of OR 95% condidence interval      | -                                            |
| `OR_95L`            | `float64`  | lower bound of OR 95% condidence interval      | -                                            |
| `HR`                | `float64`  | hazard ratio                                   | -                                            |
| `HR_95U`            | `float64`  | upper bound of HR 95% condidence interval      | -                                            |
| `HR_95L`            | `float64`  | lower bound of HR 95% condidence interval      | -                                            |
| `CHISQ`             | `float64`  | chi square                                     | -                                            |
| `Z`                 | `float64`  | z score                                        | -                                            |
| `T`                 | `float64`  | t statistics                                   | -                                            |
| `F`                 | `float64`  | F statistics                                   | -                                            |
| `P`                 | `float64`  | P value                                        | -                                            |
| `P_MANTISSA`        | `float64`  | P mantissa                                     | -                                            |
| `P_EXPONENT`        | `float64`  | P exponent                                     | -                                            |
| `MLOG10P`           | `float64`  | $-log_{10}(P)$                                 | -                                            |
| `SNPR2`             | `float64`  | per variant R2                                 | -                                            |
| `DOF`               | `Int64`    | degree of freedom                              | -                                            |
| `P_HET`             | `float64`  | heterogeneity test P value                     | -                                            |
| `I2_HET`            | `float64`  | heterogeneity I2                               | -                                            |
| `DENSITY`           | `Int64`    | signal density                                 | -                                            |
| `N`                 | `Int64`    | total sample size                              | -                                            |
| `N_CASE`            | `Int64`    | number of cases                                | -                                            |
| `N_CONTROL`         | `Int64`    | number of controls                             | -                                            |
| `GENENAME`          | `string`   | nearest gene symbol                            | -                                            |
| `CIS/TRANS`         | `string`   | whether the variant is in cis or trans region  | `Cis`,`Trans`,`NoReference`                  |
| `DISTANCE_TO_KNOWN` | `Int64`    | distance to nearest known variants             | -                                            |
| `LOCATION_OF_KNOWN` | `string`   | relative location to nearest known variants    | `Same`,`Upstream`,`Downstream`,`NoReference` |
| `KNOWN_ID`          | `string`   | nearest known variant ID                       | -                                            |
| `KNOWN_PUBMED_ID`   | `string`   | pubmed ID of the known variant                 | -                                            |
| `KNOWN_AUTHOR`      | `string`   | author of the study                            | -                                            |
| `KNOWN_SET_VARIANT` | `string`   | known set and overlapping variant              | -                                            |
| `KNOWN_VARIANT`     | `string`   | known variant overlapping with the variant     | -                                            |
| `KNOWN_SET`         | `string`   | variant set of the known variant               | -                                            |
| `NOVEL`             | `string`   | if the identified variants are novel           | -                                            |

"""