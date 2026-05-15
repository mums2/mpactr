# Mispicked ions filter

`filter_mispicked_ions()` identifies ions that were incorrectly split
into separate features during preprocessing. This filter checks the
feature table for similar ions in terms of mass and retention time.
Peaks found to be similar are merged into a single feature given
`merge_peaks` is `TRUE`.

The parameter `ringwin` is the detector saturation mass window, specific
for some instruments, such as Waters Synapse G2-Si-Q-ToF, to account for
high concentration samples.

Parameter `isowin` is the isotopic mass window, which accounts for
isotopic peaks of the same precursor mass that were incorrectly assigned
during preprocessing.

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
filter_mispicked_ions(
  mpactr_object,
  ringwin = 0.5,
  isowin = 0.01,
  trwin = 0.005,
  max_iso_shift = 3,
  merge_peaks = TRUE,
  merge_method = "sum",
  copy_object = FALSE
)
```

## Arguments

- mpactr_object:

  An `mpactr_object`. See
  [`import_data()`](https://www.mums2.org/mpactr/reference/import_data.md).

- ringwin:

  Ringing mass window or detector saturation mass window. Default = 0.5
  atomic mass units (AMU).

- isowin:

  Isotopic mass window. Default = 0.01 AMU.

- trwin:

  A `numeric` denoting the retention time threshold for assessing if
  ions should be merged. Default = 0.005.

- max_iso_shift:

  A `numeric`. Default = 3.

- merge_peaks:

  A `boolean` parameter to determine if peaks found to belong to the
  same ion should be merged in the feature table.

- merge_method:

  If merge_peaks is TRUE, a method for how similar peaks should be
  merged. Can be one of "sum".

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

data_filter <- filter_mispicked_ions(data,
  ringwin = 0.5,
  isowin = 0.01,
  trwin = 0.005,
  max_iso_shift = 3,
  merge_peaks = TRUE,
  merge_method = "sum"
)
#> ℹ Checking 1303 peaks for mispicked peaks.
#> ℹ Argument merge_peaks is: TRUE. Merging mispicked peaks with method sum.
#> ✔ 70 ions failed the mispicked filter, 1233 ions remain.
```
