# Save the Sumstats Object

GWASLab provides functions to save and load unfinished `gl.Sumstats` Objects.

## gl.dump_pickle()

Dump gl.Sumstats object.

```
gl.dump_pickle(SumstatsObject, 
                path, 
                overwrite=False)
```


## gl.load_pickle()

Load gl.Sumstats object.

```
gl.load_pickle(path)
```

## Options

| Options          | DataType      | Description                                         | Default |
|------------------|---------------|-----------------------------------------------------|---------|
| `SumstatsObject` | `gl.Sumstats` | GWASLab Sumstats Object.                            |         |
| `path`           | `string`      | File path for dumped pickle.                        |         |
| `overwrite`      | `boolean`     | If `True`, overwrite the file if it already exists. | `False` |

## Example

!!! example
    
    ```
    gl.dump_pickle(mysumstats,"./first.pickle",overwrite=True)
    
    mysumstats2 = gl.load_pickle("./first.pickle")
    ```
