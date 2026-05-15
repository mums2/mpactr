# Filter Insource ions

`filter_insource_ions()` identifies and removes in-source ion clusters
based on a Pearson correlation threshold. Groups of co-eluting features
with identical retention time are identified and used to generate
Pearson correlation matrices. Clusters with self-similarity greater than
the user-defined `cluster_threshold` within these matrices are
identified as likely belonging to a single precursor ion and is
associated insource ion. Highly correlated ions are identified and
removed.

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
filter_insource_ions(
  mpactr_object,
  cluster_threshold = 0.95,
  copy_object = FALSE
)
```

## Arguments

- mpactr_object:

  An `mpactr_object`. See
  [`import_data()`](https://mums2.github.io/mpactr/reference/import_data.md).

- cluster_threshold:

  Cluster threshold for ion deconvolution. Default = 0.95.

- copy_object:

  A `boolean` parameter that allows users to return a copied object
  instead of modifying the object.

## Value

an `mpactr_object`

## Examples

``` r
data <- import_data(example("coculture_peak_table.csv"),
  example("metadata.csv"),
  format = "Progenesis"
)

data_filter <- filter_insource_ions(data,
  cluster_threshold = 0.95
)
#> ℹ Parsing 1303 peaks for insource ions.
#> ✔ 71 ions failed the insource filter, 1232 ions remain.
```
