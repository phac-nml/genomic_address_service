# phac-nml/genomic_address_service: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.5] - 2025-03-12

### `Fixed`
- The final chunk of a file was not being processed (fixed) which resulted in query samples missing from the output [PR #25](https://github.com/phac-nml/genomic_address_service/pull/25)
- Removed max distance filtering as well because it actually results in incorrect calculations when you have multiple samples which fall outside of the threshold ranges [PR #25](https://github.com/phac-nml/genomic_address_service/pull/25)
- Fixed a bad assumption that all reference samples within a cluster would be within a range and so have added protection for missing samples [PR #25](https://github.com/phac-nml/genomic_address_service/pull/25)
- Several of these issues only occur when multiple samples were being assigned at once which was due to the fact that the lookup membership was not being updated (they had to both found new clusters and be part of the same clusters) [PR #25](https://github.com/phac-nml/genomic_address_service/pull/25)

## [0.1.4] - 2025-01-31

### `Fixed`
- Reverted change made to `assign()` class regarding cluster thresholds (> replaced with >=)
- Fixed the `assign()` class definition, `linkage_method` had been hardcoded to use only `single` even if other option was selected
- Changed some variables in the `dist_reader` and `assign` classes from static to instance to accomodate use in tests

### `Added`
- Added tests for `assign()` at different thresholds

## [0.1.3] - 2024-12-20

### `Fixed`

- Converted `data[sample_id]` to a string in the `format_df` function with `assign.py` to prevent `AttributeErrors` when non-string values are in the genomic address. [PR14](https://github.com/phac-nml/genomic_address_service/pull/14)
- Updated `buildNewick` formula to use cophenetic distances for branch lengths, aligning cluster visualization with BioNumerics dendrogram representation. [PR15](https://github.com/phac-nml/genomic_address_service/pull/15)

### `Added`

- Fixed pytests [PR7](https://github.com/phac-nml/genomic_address_service/pull/7)
- Added github actions for pytest and branch protection [PR7](https://github.com/phac-nml/genomic_address_service/pull/7)

## v1.0dev - [date]

Initial release of phac-nml/genomic_address_service

### `Added`

### `Fixed`

Changed README format to standard DAAD README, added useage arguments.

### `Dependencies`

### `Deprecated`

[0.1.3]: https://github.com/phac-nml/genomic_address_service/releases/tag/0.1.3
[0.1.4]: https://github.com/phac-nml/genomic_address_service/releases/tag/0.1.4
[0.1.5]: https://github.com/phac-nml/genomic_address_service/releases/tag/0.1.5
