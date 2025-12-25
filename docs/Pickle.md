# Save and Load the Sumstats Object

GWASLab provides functions to save and load unfinished `gl.Sumstats` objects. This is useful for:
- Saving intermediate results during long processing pipelines
- Resuming work on partially processed sumstats
- Sharing processed sumstats objects with others

## Methods

### `mysumstats.to_pickle()`

Save a Sumstats object to a pickle file. This is the recommended method as it's called directly on the Sumstats object.

```
mysumstats.to_pickle(path="~/mysumstats.pickle", 
                     overwrite=False)
```

**Parameters:**

| Parameter | DataType | Description | Default |
|-----------|------|-------------|---------|
| `path` | `str` | File path for the pickle file. Supports `~` for home directory expansion. | `"~/mysumstats.pickle"` |
| `overwrite` | `bool` | If `True`, overwrite the file if it already exists. If `False`, skip saving if file exists. | `False` |

**Returns:** None (saves the object to disk)

### `gl.dump_pickle()`

Alternative function-based approach to save a Sumstats object.

```
gl.dump_pickle(SumstatsObject, 
                path, 
                overwrite=False)
```

**Parameters:**

| Parameter | DataType | Description | Default |
|-----------|------|-------------|---------|
| `SumstatsObject` | `gl.Sumstats` | GWASLab Sumstats Object to save. | Required |
| `path` | `str` | File path for the pickle file. | Required |
| `overwrite` | `bool` | If `True`, overwrite the file if it already exists. If `False`, skip saving if file exists. | `False` |

**Returns:** None (saves the object to disk)

### `gl.load_pickle()`

Load a previously saved Sumstats object from a pickle file.

```
gl.load_pickle(path)
```

**Parameters:**

| Parameter | DataType | Description | Default |
|-----------|------|-------------|---------|
| `path` | `str` | File path to the pickle file. | Required |

**Returns:** `gl.Sumstats` object

## Examples

!!! example "Using to_pickle() method (recommended)"
    ```
    import gwaslab as gl
    
    # Load and process sumstats
    mysumstats = gl.Sumstats("sumstats.txt.gz", fmt="plink2")
    mysumstats.basic_check()
    mysumstats.fix_chr()
    mysumstats.fix_pos()
    
    # Save the processed object
    mysumstats.to_pickle("./processed_sumstats.pickle", overwrite=True)
    
    # Later, load it back
    mysumstats2 = gl.load_pickle("./processed_sumstats.pickle")
    ```

!!! example "Using dump_pickle() function"
    ```
    import gwaslab as gl
    
    # Process sumstats
    mysumstats = gl.Sumstats("sumstats.txt.gz", fmt="plink2")
    mysumstats.basic_check()
    
    # Save using function
    gl.dump_pickle(mysumstats, "./first.pickle", overwrite=True)
    
    # Load it back
    mysumstats2 = gl.load_pickle("./first.pickle")
    ```

!!! example "Save with custom path"
    ```
    # Save to home directory
    mysumstats.to_pickle("~/my_analysis/sumstats.pickle")
    
    # Save to current directory
    mysumstats.to_pickle("./sumstats.pickle", overwrite=True)
    
    # Save with absolute path
    mysumstats.to_pickle("/data/analysis/sumstats.pickle")
    ```

## Notes

- **File format**: The objects are saved using Python's `pickle` module in binary format (`.pickle` extension recommended)
- **Compatibility**: Pickled objects are generally compatible across GWASLab versions, but may require the same Python version
- **What's saved**: The entire Sumstats object is saved, including:
  - The data DataFrame
  - Metadata
  - Log history
  - All processed results (clumps, associations, etc.)
- **File existence**: By default, if a file already exists at the specified path, it will not be overwritten unless `overwrite=True`
- **Home directory**: The `~` symbol in paths is automatically expanded to your home directory
