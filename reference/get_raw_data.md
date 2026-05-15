# Return the input peak table from mpactr object.

`get_raw_data` a wrapper function to return the meta data object of the
given mpactr object.

## Usage

``` r
get_raw_data(mpactr_object)
```

## Arguments

- mpactr_object:

  The mpactr object that is created by calling the import_data()
  function.

## Value

a `data.table`.

## Examples

``` r
data <- import_data(example("coculture_peak_table.csv"),
  example("metadata.csv"),
  format = "Progenesis"
)

raw_data <- get_raw_data(data)
```
