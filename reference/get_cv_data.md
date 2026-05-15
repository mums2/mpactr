# Get CV values.

`get_cv_data()` is a wrapper function to return cv (coefficient of
variation) calculated with
[`filter_cv()`](https://www.mums2.org/mpactr/reference/filter_cv.md).

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
data <- import_data(
  example_path("coculture_peak_table.csv"),
  example_path("metadata.csv"),
  format = "Progenesis"
)

data_filter <- filter_cv(data,
  cv_threshold = 0.01,
)
#> ℹ Parsing 1303 peaks for replicability across technical replicates.
#> ✔ 1261 ions failed the cv_filter filter, 42 ions remain.

cv <- get_cv_data(data_filter)
head(cv)
#>    Compound  biological_group sample_code passes_cv_filter          cv
#>      <char>            <char>      <char>           <lgcl>       <num>
#> 1:        1            Blanks     UM1846B            FALSE  0.86736617
#> 2:        1  JC28 monoculture     UM1847B            FALSE  0.01656991
#> 3:        1   JC1 monoculture     UM1848B            FALSE -1.00000000
#> 4:        1 ANG18 monoculture     UM1849B            FALSE  0.03502637
#> 5:        1 ANGDT monoculture     UM1850B            FALSE  0.09073311
#> 6:        1         Coculture     UM1852B            FALSE  0.05760617
```
