researved_header="""

Here is a list of all headers (gl.Sumstats.data) used in GWASLab functions. These headers are fixed to avoid ambiguity within GWASLab framework.

| Header              | Keyword of args in preformat | Datatype   | Description                                    | Note                                         |
|---------------------|------------------------------|------------|------------------------------------------------|----------------------------------------------|
| `SNPID`             | snpid                        | `string`   | variant ID (CHR:POS:NEA:EA)                    | GWASLab assumes EA=ALT and NEA=REF           |
| `rsID`              | rsid                         | `string`   | dbSNP rsID                                     | -                                            |
| `CHR`               | chrom                        | `Int64`    | chromosome number (X 23, Y 24, MT 25)          | -                                            |
| `POS`               | pos                          | `Int64`    | base pair position                             | -                                            |
| `EA`                | ea                           | `category` | effect allele                                  | -                                            |
| `NEA`               | nea                          | `category` | non-effect allele                              | -                                            |
| `STATUS`            | status                       | `category` | variant standardization & harmonization status | lookup_status() to check                                            |
| `REF`               | ref                          | `category` | reference allele in reference genome           | -                                            |
| `ALT`               | alt                          | `category` | alternative allele                             | -                                            |
| `EAF`               | eaf                          | `float64`  | effect allele frequency                        | -                                            |
| `NEAF`              | neaf                         | `float64`  | non-effect allele frequency                    | -                                            |
| `MAF`               | maf                          | `float64`  | minor allele frequency                         | if EAF, use fill_data to get MAF                                            |
| `INFO`              | info                         | `float32`  | imputation INFO/RSQ                            | -                                            |
| `BETA`              | beta                         | `float64`  | effect size beta                               | -                                            |
| `SE`                | se                           | `float64`  | standard error of beta                         | -                                            |
| `BETA_95U`          | beta_95U                     | `float64`  | upper bound of beta 95% condidence interval    | -                                            |
| `BETA_95L`          | beta_95L                     | `float64`  | lower bound of beta 95% condidence interval    | -                                            |
| `OR`                | OR                           | `float64`  | odds ratio                                     | -                                            |
| `OR_95U`            | OR_95U                       | `float64`  | upper bound of OR 95% condidence interval      | -                                            |
| `OR_95L`            | OR_95L                       | `float64`  | lower bound of OR 95% condidence interval      | -                                            |
| `HR`                | HR                           | `float64`  | hazard ratio                                   | -                                            |
| `HR_95U`            | HR_95U                       | `float64`  | upper bound of HR 95% condidence interval      | -                                            |
| `HR_95L`            | HR_95L                       | `float64`  | lower bound of HR 95% condidence interval      | -                                            |
| `CHISQ`             | chisq                        | `float64`  | chi square                                     | -                                            |
| `Z`                 | z                            | `float64`  | z score                                        | -                                            |
| `T`                 | t                            | `float64`  | t statistics                                   | -                                            |
| `F`                 | f                            | `float64`  | F statistics                                   | -                                            |
| `P`                 | p                            | `float64`  | P value                                        | -                                            |
| `P_MANTISSA`        | p_mantissa                   | `float64`  | P mantissa                                     | -                                            |
| `P_EXPONENT`        | p_exponent                   | `float64`  | P exponent                                     | -                                            |
| `MLOG10P`           | mlog10p                      | `float64`  | $-log_{10}(P)$                                 | -                                            |
| `SNPR2`             | snpr2                        | `float64`  | per variant R2                                 | -                                            |
| `DOF`               | dof                          | `Int64`    | degree of freedom                              | -                                            |
| `P_HET`             | phet                         | `float64`  | heterogeneity test P value                     | -                                            |
| `I2_HET`            | i2                           | `float64`  | heterogeneity I2                               | -                                            |
| `DENSITY`           | density                      | `Int64`    | signal density                                 | -                                            |
| `N`                 | n                            | `Int64`    | total sample size                              | -                                            |
| `N_CASE`            | ncase                        | `Int64`    | number of cases                                | -                                            |
| `N_CONTROL`         | ncontrol                     | `Int64`    | number of controls                             | -                                            |
| `GENENAME`          |                       | `string`   | nearest gene symbol                            | -                                            |
| `CIS/TRANS`         |                      | `string`   | whether the variant is in cis or trans region  | `Cis`,`Trans`,`NoReference`                  |
| `DISTANCE_TO_KNOWN` |              | `Int64`    | distance to nearest known variants             | -                                            |
| `LOCATION_OF_KNOWN` |              | `string`   | relative location to nearest known variants    | `Same`,`Upstream`,`Downstream`,`NoReference` |
| `KNOWN_ID`          |                       | `string`   | nearest known variant ID                       | -                                            |
| `KNOWN_PUBMED_ID`   |                | `string`   | pubmed ID of the known variant                 | -                                            |
| `KNOWN_AUTHOR`      |                   | `string`   | author of the study                            | -                                            |
| `KNOWN_SET_VARIANT` |              | `string`   | known set and overlapping variant              | -                                            |
| `KNOWN_VARIANT`     |                  | `string`   | known variant overlapping with the variant     | -                                            |
| `KNOWN_SET`         |                      | `string`   | variant set of the known variant               | -                                            |
| `NOVEL`             |                          | `string`   | if the identified variants are novel           | -                                            |

"""

researved_header_json="""
{
  "SNPID":{"Keyword":"snpid","Datatype":"string","Description":"variant ID (CHR:POS:NEA:EA)","Note":"GWASLab assumes EA=ALT and NEA=REF"},
  "rsID":{"Keyword":"rsid","Datatype":"string","Description":"dbSNP rsID","Note":"-"},
  "CHR":{"Keyword":"chrom","Datatype":"Int64","Description":"chromosome number (X 23, Y 24, MT 25)","Note":"-"},
  "POS":{"Keyword":"pos","Datatype":"Int64","Description":"base pair position","Note":"-"},
  "EA":{"Keyword":"ea","Datatype":"category","Description":"effect allele","Note":"-"},
  "NEA":{"Keyword":"nea","Datatype":"category","Description":"non-effect allele","Note":"-"},
  "STATUS":{"Keyword":"status","Datatype":"category","Description":"7-digit variant status","Note":"lookup_status() to check"},
  "REF":{"Keyword":"ref","Datatype":"category","Description":"reference allele in reference genome","Note":"-"},
  "ALT":{"Keyword":"alt","Datatype":"category","Description":"alternative allele","Note":"-"},
  "EAF":{"Keyword":"eaf","Datatype":"float64","Description":"effect allele frequency","Note":"-"},
  "NEAF":{"Keyword":"neaf","Datatype":"float64","Description":"non-effect allele frequency","Note":"-"},
  "MAF":{"Keyword":"maf","Datatype":"float64","Description":"minor allele frequency","Note":"if EAF, use fill_data to get MAF"},
  "INFO":{"Keyword":"info","Datatype":"float32","Description":"imputation INFO/RSQ","Note":"-"},
  "BETA":{"Keyword":"beta","Datatype":"float64","Description":"effect size beta","Note":"-"},
  "SE":{"Keyword":"se","Datatype":"float64","Description":"standard error of beta","Note":"-"},
  "BETA_95U":{"Keyword":"beta_95U","Datatype":"float64","Description":"upper bound of beta 95% condidence interval","Note":"-"},
  "BETA_95L":{"Keyword":"beta_95L","Datatype":"float64","Description":"lower bound of beta 95% condidence interval","Note":"-"},
  "OR":{"Keyword":"OR","Datatype":"float64","Description":"odds ratio","Note":"-"},
  "OR_95U":{"Keyword":"OR_95U","Datatype":"float64","Description":"upper bound of OR 95% condidence interval","Note":"-"},
  "OR_95L":{"Keyword":"OR_95L","Datatype":"float64","Description":"lower bound of OR 95% condidence interval","Note":"-"},
  "HR":{"Keyword":"HR","Datatype":"float64","Description":"hazard ratio","Note":"-"},
  "HR_95U":{"Keyword":"HR_95U","Datatype":"float64","Description":"upper bound of HR 95% condidence interval","Note":"-"},
  "HR_95L":{"Keyword":"HR_95L","Datatype":"float64","Description":"lower bound of HR 95% condidence interval","Note":"-"},
  "CHISQ":{"Keyword":"chisq","Datatype":"float64","Description":"chi square","Note":"-"},
  "Z":{"Keyword":"z","Datatype":"float64","Description":"z score","Note":"-"},
  "T":{"Keyword":"t","Datatype":"float64","Description":"t statistics","Note":"-"},
  "F":{"Keyword":"f","Datatype":"float64","Description":"F statistics","Note":"-"},
  "P":{"Keyword":"p","Datatype":"float64","Description":"P value","Note":"-"},
  "P_MANTISSA":{"Keyword":"p_mantissa","Datatype":"float64","Description":"P mantissa","Note":"-"},
  "P_EXPONENT":{"Keyword":"p_exponent","Datatype":"float64","Description":"P exponent","Note":"-"},
  "MLOG10P":{"Keyword":"mlog10p","Datatype":"float64","Description":"-log10(P)","Note":"-"},
  "SNPR2":{"Keyword":"snpr2","Datatype":"float64","Description":"per variant R2","Note":"-"},
  "DOF":{"Keyword":"dof","Datatype":"Int64","Description":"degree of freedom","Note":"-"},
  "P_HET":{"Keyword":"phet","Datatype":"float64","Description":"heterogeneity test P value","Note":"-"},
  "I2_HET":{"Keyword":"i2","Datatype":"float64","Description":"heterogeneity I2","Note":"-"},
  "DENSITY":{"Keyword":"density","Datatype":"Int64","Description":"signal density","Note":"-"},
  "N":{"Keyword":"n","Datatype":"Int64","Description":"total sample size","Note":"-"},
  "N_CASE":{"Keyword":"ncase","Datatype":"Int64","Description":"number of cases","Note":"-"},
  "N_CONTROL":{"Keyword":"ncontrol","Datatype":"Int64","Description":"number of controls","Note":"-"},
  "GENENAME":{"Keyword":"","Datatype":"string","Description":"nearest gene symbol","Note":"-"},
  "CIS/TRANS":{"Keyword":"","Datatype":"string","Description":"whether the variant is in cis or trans region","Note":"Cis,Trans,NoReference"},
  "DISTANCE_TO_KNOWN":{"Keyword":"","Datatype":"Int64","Description":"distance to nearest known variants","Note":"-"},
  "LOCATION_OF_KNOWN":{"Keyword":"","Datatype":"string","Description":"relative location to nearest known variants","Note":"Same,Upstream,Downstream,NoReference"},
  "KNOWN_ID":{"Keyword":"","Datatype":"string","Description":"nearest known variant ID","Note":"-"},
  "KNOWN_PUBMED_ID":{"Keyword":"","Datatype":"string","Description":"pubmed ID of the known variant","Note":"-"},
  "KNOWN_AUTHOR":{"Keyword":"","Datatype":"string","Description":"author of the study","Note":"-"},
  "KNOWN_SET_VARIANT":{"Keyword":"","Datatype":"string","Description":"known set and overlapping variant","Note":"-"},
  "KNOWN_VARIANT":{"Keyword":"","Datatype":"string","Description":"known variant overlapping with the variant","Note":"-"},
  "KNOWN_SET":{"Keyword":"","Datatype":"string","Description":"variant set of the known variant","Note":"-"},
  "NOVEL":{"Keyword":"","Datatype":"string","Description":"if the identified variants are novel","Note":"-"}
}
"""
import json
def _build_reserved_header(j):
    d = json.loads(j)
    rows = []
    rows.append("Here is a list of all headers (gl.Sumstats.data) used in GWASLab functions. These headers are fixed to avoid ambiguity within GWASLab framework.")
    rows.append("")
    rows.append("| Header              | Keyword of args in preformat | Datatype   | Description                                    | Note                                         |")
    rows.append("|---------------------|------------------------------|------------|------------------------------------------------|----------------------------------------------|")
    for k, v in d.items():
        keyword = v.get("Keyword","")
        datatype = v.get("Datatype","")
        description = v.get("Description","")
        note = v.get("Note","")
        if keyword == "":
            kw = " "
        else:
            kw = keyword
        if datatype == "":
            dt = " "
        else:
            dt = f"`{datatype}`"
        rows.append(f"| `{k}` | {kw} | {dt} | {description} | {note} |")
    rows.append("")
    return "\n".join(rows)
researved_header = _build_reserved_header(researved_header_json)
