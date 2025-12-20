# Liftover

GWASLab can directly liftover the positions of variants in sumstats.

## .liftover()

```
mysumstats.liftover(threads=3, 
                    from_build="19", 
                    to_build="38",
                    remove=True)
```

!!! note
    GWASLab will only liftover basepair positions. If needed, please perform harmonization using the reference file of the target build.  

## Options

| `.liftover()` options | DataType         | Description                            | Default |
|-----------------------|------------------|----------------------------------------|---------|
| `threads`             | `int`            | Number of threads to use for liftover. | `1`     |
| `from_build`          | `"19"` or `"38"` | Original genome build                  | -       |
| `to_build`            | `"19"` or `"38"` | Target genome build                    | -       |
| `remove`              | `boolean`        | If True, remove unmapped variants      | `True`  |

!!! quote
    This method is based on [GitHub - jeremymcrae/liftover: liftover for python, made fast with cython](https://github.com/jeremymcrae/liftover)

## Example

!!!example
    
    ```
    mysumstats.liftover(threads=3, from_build="19", to_build="38")
    
    ```
    

