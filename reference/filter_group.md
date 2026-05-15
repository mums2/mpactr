# Filter Ions by Group

Filter Ions by Group

## Usage

``` r
filter_group(
  mpactr_object,
  group_threshold = 0.01,
  group_to_remove,
  remove_ions = TRUE,
  copy_object = FALSE
)
```

## Arguments

- mpactr_object:

  An `mpactr_object`. See
  [`import_data()`](https://mums2.github.io/mpactr/reference/import_data.md).

- group_threshold:

  Relative abundance threshold at which to remove ions. Default = 0.01.

- group_to_remove:

  Biological group name to remove ions from.

- remove_ions:

  A `boolean` parameter. If `TRUE` failing ions will be removed from the
  peak table. Default = TRUE.

- copy_object:

  A `boolean` parameter that allows users to return a copied object
  instead of modifying the object.

## Value

an `mpactr_object`.

## Details

`filter_group()` removes feature ions that are present in a user-defined
group based on a relative abundance threshold. This could be
particularly useful to filter out features found present in solvent
blank samples. Further, this filter can be ultilized to remove features
in media blank sample for experiments on microbial cultures. The
presence or absence of features in a group of samples is determined by
first averaging injection replicates and then averaging biological
replicates within each biological treatment group. A feature is present
in a group if its abundance is greater than the user-defined
`group_threshold`. The default is 0.01, meaning a feature is removed if
its abundance is 1% of that in the sample group in which it is most
abundant. For example, blank filtering can remove features whose mean
abundance in solvent blank injections is greater than 1% of their
maximum mean abundance in experimental samples.

If you would like to remove features found in media blank samples, we
recommend testing the `group_threshold` parameter.

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

## Examples

``` r
data <- import_data(example("coculture_peak_table.csv"),
  example("metadata.csv"),
  format = "Progenesis"
)

data_filter <- filter_group(data,
  group_threshold = 0.01,
  group_to_remove = "Blanks",
  remove_ions = TRUE
)
#> ℹ Parsing 1303 peaks based on the sample group: Blanks.
#> ℹ Argument remove_ions is: TRUE.Removing peaks from Blanks.
#> ✔ 796 ions failed the Blanks filter, 507 ions remain.
```
