# phac-nml/genomic_address_service: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.3] - 2024/11/..

### `Fixed`

- Converted `data[sample_id]` to a string in the `format_df` function with `assign.py` to prevent `AttributeErrors` when non-string values are in the genomic address.

## v1.0dev - [date]

Initial release of phac-nml/genomic_address_service

### `Added`

### `Fixed`

Changed README format to standard DAAD README, added useage arguments.

### `Dependencies`

### `Deprecated`

[0.1.3]: https://github.com/phac-nml/genomic_address_service/releases/tag/0.1.3
