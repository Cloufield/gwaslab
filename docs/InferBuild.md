# Infer Genome Build

GWASLab uses chromosome and base pair position information from built-in Hapmap3 SNPs to infer the reference genome build (GRCh37/hg19 or GRCh38/hg38) for your summary statistics.

## How It Works

The function compares the CHR:POS coordinates in your sumstats against Hapmap3 SNP coordinates for both hg19 and hg38. The genome build with the most matching variants is inferred as the correct build version.

1. **Loads Hapmap3 reference data** for both hg19 and hg38
2. **Matches CHR:POS coordinates** between your sumstats and reference data
3. **Counts matches** for each build version
4. **Infers the build** based on which has more matches
5. **Updates STATUS codes** to reflect the inferred build (if `change_status=True`)
6. **Updates metadata** in the Sumstats object

## .infer_build()

```
mysumstats.infer_build(verbose=True)
```

## Parameters

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `verbose` | `boolean` | If True, print detailed logs including match counts for each build | `True` |

!!! note "Reliability"
    The results are more reliable if you have a reasonable amount of variants. The function will issue a warning if fewer than 10,000 matching variants are found, as this may indicate unreliable inference.

!!! note "STATUS Code Updates"
    The function automatically updates the STATUS code to reflect the inferred build:
    - For hg19: Sets digit 1 to "1" (build = 19)
    - For hg38: Sets digit 1 to "3" and digit 2 to "8" (build = 38)
    
    This ensures the STATUS code accurately reflects the inferred genome build.

## Return Value

The function updates the Sumstats object in place:
- Updates `mysumstats.data` with modified STATUS codes
- Updates `mysumstats.build` with the inferred build ("19" or "38")
- Updates `mysumstats.meta["gwaslab"]["genome_build"]` with the inferred build

## Examples

!!! example "Basic usage"
    ```
    # Infer genome build
    mysumstats.infer_build()
    
    # Check the inferred build
    print(mysumstats.build)  # "19" or "38"
    ```


!!! example "Quiet mode"
    ```
    # Infer build without verbose output
    mysumstats.infer_build(verbose=False)
    ```

## Output

The function provides detailed logging output:

```
 -Loading Hapmap3 variants data...
 -CHR and POS will be used for matching...
 -Matching variants for hg19: num_hg19 = 12345
 -Matching variants for hg38: 5678
 -Since num_hg19 >> num_hg38, set the genome build to hg19 for the STATUS code....
```

If fewer than 10,000 matches are found:
```
⚠️ Please be cautious due to the limited number of variants.
```

## Notes

- **Requires CHR and POS columns**: The function requires valid chromosome and position information in your sumstats
- **Hapmap3 reference**: Uses built-in Hapmap3 SNP coordinates for matching
- **Automatic metadata update**: The inferred build is automatically stored in the Sumstats object's metadata
- **STATUS code integration**: The STATUS code is automatically updated to reflect the inferred build, which is important for downstream harmonization steps
- **Uncertain inference**: If the match counts are equal, the function cannot determine the build and will log a message

## When to Use

- When the genome build of your sumstats is unknown
- Before performing harmonization or liftover operations
- To verify the genome build of your data
- As part of quality control workflows

## Related Functions

- `.set_build()`: Manually set the genome build
- `.liftover()`: Convert coordinates between genome builds (requires knowing the source build)
