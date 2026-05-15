# Get CV values.

`get_cv_data()` is a wrapper function to return cv (coefficient of
variation) calculated with
[`filter_cv()`](https://mums2.github.io/mpactr/reference/filter_cv.md).

## Usage

``` r
get_cv_data(mpactr_object)
```

## Arguments

- mpactr_object:

  The mpactr object that is created by calling the import_data()
  function.

## Value

a `data.table` reporting the mean and median coefficient of variation
for each input ion.

## Examples

``` r
data <- import_data(example("coculture_peak_table.csv"),
  example("metadata.csv"),
  format = "Progenesis"
)

data_filter <- filter_cv(data,
  cv_threshold = 0.01,
  cv_param = "median"
)
#> ℹ Parsing 1303 peaks for replicability across technical replicates.
#> ✔ 1298 ions failed the cv_filter filter, 5 ions remain.

cv <- get_cv_data(data_filter)
head(cv)
#>    Compound    mean_cv  median_cv
#>      <char>      <num>      <num>
#> 1:        1 0.21346035 0.05760617
#> 2:        2 0.91799997 0.91799997
#> 3:        3 0.90848851 0.90848851
#> 4:        4 0.06782583 0.06782583
#> 5:        5 0.02727037 0.02727037
#> 6:        6 0.09005716 0.09005716
```
