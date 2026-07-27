# Changelog


---

## [Unreleased]
### Changed
- **Renamed the package from `pairs_to_parquet` to `pairtools_parquet`**, along
  with the CLI entry point and the `@PG` provenance records it writes. This is a
  breaking change for imports and for the command name.

### Fixed
- `pairtools_parquet.lib` was missing an `__init__.py` and so was dropped from
  non-editable installs.
- Declared the `pandas` and `numpy` runtime dependencies; dropped the unused
  `fire` dependency.
- Tests no longer write into tracked files under `tests/data/`.

---

## [0.2.0] - 2025-11-26
### Added
- Includes CLI tool for selecting `.parquet` files only.
- Added testing.

---

## [0.1.0] - 2025-10-21
### Added
- First prototype version for internal testing of **pairs-to-parquet**.
- Includes CLI tool for sorting `.pairs` and `.parquet` files.
- Added documentation and optional dependencies for testing and docs.

