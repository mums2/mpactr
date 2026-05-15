# Changelog

## mpactr (development version)

## mpactr 0.3.2

CRAN release: 2026-05-14

- Updated mpactr so that all column names inside the metadata file are
  lowercase. You may still enter you data with uppercase column names,
  but we will force all data inside the metadata file to be lowercase.
  This change happens within the
  [`import_data()`](https://www.mums2.org/mpactr/reference/import_data.md)
  function.
- Updated vignette to reflect current changes.
- Renamed the `get_meta_data()` to
  [`get_metadata()`](https://www.mums2.org/mpactr/reference/get_metadata.md).

## mpactr 0.3.1

CRAN release: 2025-09-27

- Fixed the Valgrind memory issues. This package should now work with
  valgrind, clang-asan, and gcc-asan.

- Updated documentation for
  [`filter_cv()`](https://www.mums2.org/mpactr/reference/filter_cv.md)
  to reflect current changes.

- Updated vignette to reflect
  [`filter_cv()`](https://www.mums2.org/mpactr/reference/filter_cv.md)
  changes.

- Added a
  [`limit_cores()`](https://www.mums2.org/mpactr/reference/limit_cores.md)
  function to limit cores used when going through CRAN checks.

## mpactr 0.3.0

CRAN release: 2025-09-07

- Major changes to
  [`filter_cv()`](https://www.mums2.org/mpactr/reference/filter_cv.md).
  Old
  [`filter_cv()`](https://www.mums2.org/mpactr/reference/filter_cv.md)
  would calculate the coefficient of variance of the biological groups,
  calculate the average/median, and remove the feature if the
  average/mean was above the threshold. To account for some groups
  having high variation by default, the new
  [`filter_cv()`](https://www.mums2.org/mpactr/reference/filter_cv.md)
  will now zero out groups if they are above the threshold. If every
  group in the feature has been zeroed out, we will remove the feature.

- [`filter_cv()`](https://www.mums2.org/mpactr/reference/filter_cv.md)
  has been optimized to use c++ code, so there should be significant
  speed improvements.

- You are also now able to import `data.frames` for peak table input in
  [`import_data()`](https://www.mums2.org/mpactr/reference/import_data.md).

## mpactr 0.2.1

CRAN release: 2025-03-28

## mpactr 0.1.0

CRAN release: 2024-09-10

## mpactr 0.1.0

CRAN release: 2024-09-10

- Initial CRAN submission.

- Contains four main filters:
  [`filter_mispicked_ions()`](https://www.mums2.org/mpactr/reference/filter_mispicked_ions.md),
  [`filter_group()`](https://www.mums2.org/mpactr/reference/filter_group.md),
  [`filter_cv()`](https://www.mums2.org/mpactr/reference/filter_cv.md),
  and
  [`filter_insource_ions()`](https://www.mums2.org/mpactr/reference/filter_insource_ions.md)
  for the correction of tandem MS/MS peaks.
