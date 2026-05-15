# Summary of Filtering

Parses an mpactr object and extracts a summary of all applied filters.
Specifically, the fate of each input ion is reported as ion status.
Status options are: Passed, mispicked, group, replicability, and
insouce. A status of Passed ions is returned for ions that passed all
applied filters and therefore are expected to be high quality ions. Ions
tagged as group, mispicked, replicability, or ionsource were removed
during the corresponding filter.

## Usage

``` r
qc_summary(mpactr_object)
```

## Arguments

- mpactr_object:

  an `mpactr_object`.

## Value

a `data.table` reporting the number of high quality ions ("Passed") or
the filter in which they were removed.

## Examples

``` r
data <- import_data(
  example_path("coculture_peak_table.csv"),
  example_path("metadata.csv"),
  format = "Progenesis"
)

data_filter <- filter_mispicked_ions(data,
  ringwin = 0.5,
  isowin = 0.01,
  trwin = 0.005,
  max_iso_shift = 3,
  merge_peaks = TRUE
)
#> ℹ Checking 1303 peaks for mispicked peaks.
#> ℹ Argument merge_peaks is: TRUE. Merging mispicked peaks with method sum.
#> ✔ 70 ions failed the mispicked filter, 1233 ions remain.

summary <- qc_summary(data_filter)
summary
#>          status compounds
#>          <char>    <char>
#>    1:    Passed         1
#>    2:    Passed        10
#>    3:    Passed       100
#>    4:    Passed      1000
#>    5:    Passed      1001
#>   ---                    
#> 1299: mispicked      1267
#> 1300: mispicked       546
#> 1301: mispicked       495
#> 1302: mispicked      1289
#> 1303: mispicked      1290
```
