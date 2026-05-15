# Filter Non-reproducible ions

`filter_cv()` removes feature ions that are found to be non-reproducible
between technical injection replicates. Reproducibility is assessed via
mean or median coefficient of variation (CV) between technical
replicates. As such, this filter is expecting an input dataset with at
least two replicate injections per sample.

`copy_object`: mpactr is built on an R6 class-system, meaning it
operates on reference semantics in which data is updated *in-place*.
Compared to a shallow copy, where only data pointers are copied, or a
deep copy, where the entire data object is copied in memory, any changes
to the original data object, regardless if they are assigned to a new
object, result in changes to the original data object. We recommend
using the default `copy_object = FALSE` as this makes for an extremely
fast and memory-efficient way to chain mpactr filters together; however,
if you would like to run the filters individually with traditional R
style objects, you can set `copy_object` to `TRUE` as shown in the
filter examples.

## Usage

``` r
filter_cv(mpactr_object, cv_threshold = NULL, cv_param, copy_object = FALSE)
```

## Arguments

- mpactr_object:

  An `mpactr_object`. See
  [`import_data()`](https://mums2.github.io/mpactr/reference/import_data.md).

- cv_threshold:

  Coefficient of variation threshold. A lower cv_threshold will result
  in more stringent filtering and higher reproducibility. Recommended
  values between 0.2 - 0.5.

- cv_param:

  Coefficient of variation (CV) statistic to use for filtering Options
  are "mean" or "median", corresponding to mean and median CV,
  respectively.

- copy_object:

  A `boolean` parameter that allows users to return a copied object
  instead of modifying the object.

## Value

an `mpactr_object`.

## Examples

``` r

data <- import_data(example("coculture_peak_table.csv"),
  example("metadata.csv"),
  format = "Progenesis"
)

data_filter <- filter_cv(data,
  cv_threshold = 0.01,
  cv_param = "mean",
  copy_object = TRUE
)
#> ℹ Parsing 1303 peaks for replicability across technical replicates.
#> ✔ 1298 ions failed the cv_filter filter, 5 ions remain.

data_filter <- filter_cv(data,
  cv_threshold = 0.01,
  cv_param = "median",
  copy_object = TRUE
)
#> ℹ Parsing 1303 peaks for replicability across technical replicates.
#> ✔ 1298 ions failed the cv_filter filter, 5 ions remain.
```
