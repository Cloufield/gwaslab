# Liftover

GWASLab can directly liftover the positions of variants in sumstats using a fast chain-based implementation.

## .liftover()

```
mysumstats.liftover(from_build="19", 
                    to_build="38",
                    remove=True)
```

!!! note
    GWASLab will only liftover basepair positions. If needed, please perform harmonization using the reference file of the target build.  

## Implementation Details

The liftover implementation uses the `sumstats_liftover` package, which directly parses UCSC chain files and performs vectorized coordinate conversion for optimal performance. 

**Key features:**
- **Fast chain-based conversion**: Directly parses UCSC chain files with vectorized operations
- **Built-in chain files**: Automatically uses built-in chain files for hg19↔hg38 conversions (downloaded from UCSC if not available)
- **Chromosome mismatch detection**: Variants mapped to different chromosomes are automatically detected and treated as unmapped

## Options

| `.liftover()` options | DataType         | Description                            | Default |
|-----------------------|------------------|----------------------------------------|---------|
| `from_build`          | `"19"` or `"38"` | Original genome build                  | -       |
| `to_build`            | `"19"` or `"38"` | Target genome build                    | -       |
| `chain_path`          | `str`            | Path to UCSC chain file (optional). If provided, `from_build` and `to_build` are optional | `None`  |
| `remove`              | `boolean`        | If True, remove unmapped variants (including chromosome mismatches) | `True`  |

!!! quote
    This method uses the [sumstats_liftover](https://github.com/statgen/sumstats_liftover) package, which provides fast chain-based coordinate conversion using vectorized operations.

## Example

!!!example
    
    ```
    # Basic usage with build numbers
    mysumstats.liftover(from_build="19", to_build="38")
    
    # Using a custom chain file
    mysumstats.liftover(chain_path="/path/to/chain/file.chain")
    
    # Keep unmapped variants in the output
    mysumstats.liftover(from_build="19", to_build="38", remove=False)
    ```
    

