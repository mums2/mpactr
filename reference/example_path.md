# Get file paths for examples

mpactr contains a number of example files in the `inst/extdata`
directory. This function makes them accessible in documentation that
shows how file paths are used in function examples.

## Usage

``` r
example_path(file = NULL)
```

## Arguments

- file:

  Name of a file. If `NULL`, all examples files will be listed.

## Value

A file path to example data stored in the `inst/extdata` directory of
the package.

## Examples

``` r
example_path()
#>  [1] "PTY087I2_dataset.csv"               "PTY087I2_extractmetadata.csv"      
#>  [3] "PTY087I2_samplelist.csv"            "coculture_peak_table.csv"          
#>  [5] "cultures_metaboscape_metadata.csv"  "cultures_metaboscape_peaktable.csv"
#>  [7] "cultures_metadata.csv"              "cultures_peak_table.csv"           
#>  [9] "ion_masses"                         "metadata.csv"                      

example_path("metadata.csv")
#> [1] "/home/runner/work/_temp/Library/mpactr/extdata/metadata.csv"
```
