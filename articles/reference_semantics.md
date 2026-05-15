# Reference Semantics

``` r

library(mpactr)
#> Loading required package: data.table
#> 
#> Attaching package: 'data.table'
#> The following object is masked from 'package:base':
#> 
#>     %notin%
#> 
#> Attaching package: 'mpactr'
#> The following object is masked from 'package:utils':
#> 
#>     example
```

## Reference semantics

mpactr is built on an R6 class-system, meaning it operates on reference
semantics in which data is updated *in-place*. Compared to a shallow
copy, where only data pointers are copied, or a deep copy, where the
entire data object is copied in memory, any changes to the original data
object, regardless if they are assigned to a new object, result in
changes to the original data object. We can see this below.

``` r

data2 <- import_data(example("cultures_peak_table.csv"),
  example("cultures_metadata.csv"),
  format = "Progenesis"
)


get_peak_table(data2)[, 1:5]
#>       Compound       mz         rt 102623_UM1848B_JC1_69_1_5004
#>         <char>    <num>      <num>                        <num>
#>    1:        1 256.0883  0.7748333                         0.00
#>    2:        2 484.2921  0.7756667                       546.56
#>    3:        3 445.2276  0.7786667                         0.00
#>    4:        4 354.1842  0.7786667                         0.00
#>    5:        5 353.1995  0.7816667                         0.00
#>   ---                                                          
#> 1330:     1330 538.3182 11.4396667                         0.00
#> 1331:     1331 424.2710 11.4315000                         0.00
#> 1332:     1332 422.1770 11.4568333                     10655.27
#> 1333:     1333 422.1776 11.4528333                       923.77
#> 1334:     1334 538.3981 11.4811667                      1176.19
#>       102623_UM1846B_Media_67_1_5005
#>                                <num>
#>    1:                           0.00
#>    2:                       16389.28
#>    3:                       22515.28
#>    4:                        6086.35
#>    5:                        5923.96
#>   ---                               
#> 1330:                           0.00
#> 1331:                           0.00
#> 1332:                        5737.01
#> 1333:                           0.00
#> 1334:                        1353.96
```

Where the raw data object has 1334 ions in the feature table.

``` r

data2_mispicked <- filter_mispicked_ions(data2,
  ringwin = 0.5,
  isowin = 0.01, trwin = 0.005,
  max_iso_shift = 3, merge_peaks = TRUE,
  merge_method = "sum",
  copy_object = FALSE
)
#> ℹ Checking 1334 peaks for mispicked peaks.
#> ℹ Argument merge_peaks is: TRUE. Merging mispicked peaks with method sum.
#> ✔ 73 ions failed the mispicked filter, 1261 ions remain.

get_peak_table(data2_mispicked)[, 1:5]
#> Key: <Compound, mz, kmd, rt>
#>       Compound       mz     kmd        rt 102423_Blank_77_1_5095
#>         <char>    <num>   <num>     <num>                  <num>
#>    1:        1 256.0883 0.08831 0.7748333                      0
#>    2:       10 340.2040 0.20399 0.7916667                      0
#>    3:      100 557.1519 0.15191 3.6925000                      0
#>    4:     1000 278.0638 0.06382 5.5228333                      0
#>    5:     1001 296.0736 0.07365 5.5246667                      0
#>   ---                                                           
#> 1257:      995 561.2726 0.27255 5.4810000                      0
#> 1258:      996 228.1430 0.14305 5.4818333                      0
#> 1259:      997 425.1873 0.18726 5.4640000                      0
#> 1260:      998 337.1987 0.19873 5.4818333                      0
#> 1261:      999 640.3299 0.32993 5.4596667                      0
```

Running the `filter_mispicked_ions` filter, with default setting
`copy_object = FALSE` (operates on reference semantics) results in 1261
ions in the feature table.

Even though we created an object called `data2_mispicked`, the original
`data2` object was also updated and now has 1261 ions in the feature
table:

``` r

get_peak_table(data2)[, 1:5]
#> Key: <Compound, mz, kmd, rt>
#>       Compound       mz     kmd        rt 102423_Blank_77_1_5095
#>         <char>    <num>   <num>     <num>                  <num>
#>    1:        1 256.0883 0.08831 0.7748333                      0
#>    2:       10 340.2040 0.20399 0.7916667                      0
#>    3:      100 557.1519 0.15191 3.6925000                      0
#>    4:     1000 278.0638 0.06382 5.5228333                      0
#>    5:     1001 296.0736 0.07365 5.5246667                      0
#>   ---                                                           
#> 1257:      995 561.2726 0.27255 5.4810000                      0
#> 1258:      996 228.1430 0.14305 5.4818333                      0
#> 1259:      997 425.1873 0.18726 5.4640000                      0
#> 1260:      998 337.1987 0.19873 5.4818333                      0
#> 1261:      999 640.3299 0.32993 5.4596667                      0
```

We recommend using the default `copy_object = FALSE` as this makes for
an extremely fast and memory-efficient way to chain mpactr filters
together (see the Filter article); however, if you would like to run the
filters individually with traditional R style objects, you can set
`copy_object` to `TRUE` as shown in the filter examples.

## Deep copy: specifics of the R6 class

The R6 class-system operates on reference semantics in which data is
updated *in-place*. Compared to a shallow copy, where only data pointers
are copied, or a deep copy, where the entire data object is copied in
memory, any changes to the original data object, regardless if they are
assigned to a new object, result in changes to the original data object.
R6 Accomplishes this by taking advantaging of the environment system in
R. Inside R, everything is created inside a base R environment. This
contains all functions, saved variables, libraries, references, etc.
Using R6 classes allows us to easily add this functionality to our R
package.

In general, R relies on reference semantics to store data away from the
outside because R environments are a container for a copious amount of
data. In a normal R session, the base R environment is the outermost
environment, allowing you to access to everything you need.

Reference semantics become noticeable when you send an environmental
variable to a function. In R, functions rely on call-by-value semantics.
Call-by-value is described as functions treating parameterized values
(values specified when calling the function) as local variables when in
the function. Anything you do to the variable in the function will have
no effect on the variables outside. This follows traditional copy by
value semantics. However, R does not allow you to send over variables by
reference due to this idea. So, you can think of functions as temporary
environments in R. What makes these environments so powerful, is the
fact that you can send an environment to a function, and it will not
copy the environment. This allows you to send variables by reference to
functions. R6 classes rely on this, and mpactr uses this for speedy
execution.

Memory usage really shines when you use R6 classes vs. a traditional
workflow, such as copy by value. In a traditional workflow, all of the
data must be copied to call functions and compute operations, using R6
classes we can minimize that problem, improving performace for large
datasets.
