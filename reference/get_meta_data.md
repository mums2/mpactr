# Return the meta_data data.table from the mpactr object.

`get_meta_data()` a wrapper function to return the meta data object of
the given mpactr object.

## Usage

``` r
get_meta_data(mpactr_object)
```

## Arguments

- mpactr_object:

  The mpactr object that is created by calling the import_data()
  function.

## Value

a `data.table`.

## Examples

``` r
data <- import_data(
  example_path("coculture_peak_table.csv"),
  example_path("metadata.csv"),
  format = "Progenesis"
)

meta_data <- get_meta_data(data)
```
