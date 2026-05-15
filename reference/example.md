# Get file paths for examples

mpactr contains a number of example files in the `inst/extdata`
directory. This function makes them accessible in documentation that
shows how file paths are used in function examples.

## Usage

``` r
example(file = NULL)
```

## Arguments

- file:

  Name of a file. If `NULL`, all examples files will be listed.

## Value

A file path to example data stored in the `inst/extdata` directory of
the package.

## Examples

``` r
example()
#> [1] "coculture_peak_table.csv"           "cultures_metaboscape_metadata.csv" 
#> [3] "cultures_metaboscape_peaktable.csv" "cultures_metadata.csv"             
#> [5] "cultures_peak_table.csv"            "ion_masses"                        
#> [7] "metadata.csv"                      

example("metadata.csv")
#> [1] "/home/runner/work/_temp/Library/mpactr/extdata/metadata.csv"
```
