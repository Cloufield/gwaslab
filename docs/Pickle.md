# Save the Sumstats Object

GWASLab provides functions to save and load unfinished `gl.Sumstats` Objects.

## Usage

```
gl.dump_pickle(SumstatsObject, path, overwrite=False)

gl.load_pickle(path)
```

## Options
- `SumstatsObject` : `gl.Sumstats()`. GWASLab Sumstats Object.
- `path` : `string`. dumped pickle file path.
- `overwrite` : `boolean`. If True, overwrite the file if it already exists.

## Example

!!! example

```
gl.dump_pickle(mysumstats,"./first.pickle",overwrite=True)
Wed Jan 11 23:29:47 2023 Start to dump the Sumstats Object.
Wed Jan 11 23:29:47 2023  -Dump the Sumstats Object to :  ./first.pickle

my2sumstats = gl.load_pickle("./first.pickle")
Wed Jan 11 23:45:59 2023 Loaded dumped Sumstats object from :  ./first.pickle
```