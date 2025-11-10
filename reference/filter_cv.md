# Filter Non-reproducible ions

`filter_cv()` removes feature ions that are found to be non-reproducible
between technical injection replicates. Reproducibility is assessed via
mean coefficient of variation (CV) between technical replicates. As
such, this filter is expecting an input dataset with at least two
replicate injections per sample.

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
filter_cv(mpactr_object, cv_threshold = NULL, copy_object = FALSE)
```

## Arguments

- mpactr_object:

  An `mpactr_object`. See
  [`import_data()`](https://www.mums2.org/mpactr/reference/import_data.md).

- cv_threshold:

  Coefficient of variation threshold. A lower cv_threshold will result
  in more stringent filtering and higher reproducibility. Recommended
  values between 0.2 - 0.5.

- copy_object:

  A `boolean` parameter that allows users to return a copied object
  instead of modifying the object.

## Value

an `mpactr_object`.

## Examples

``` r
data <- import_data(
  example_path("coculture_peak_table.csv"),
  example_path("metadata.csv"),
  format = "Progenesis"
)

data_filter <- filter_cv(data,
  cv_threshold = 0.01,
  copy_object = TRUE
)
#> ℹ Parsing 1303 peaks for replicability across technical replicates.
#> ✔ 1261 ions failed the cv_filter filter, 42 ions remain.
```
