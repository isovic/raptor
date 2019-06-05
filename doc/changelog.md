# Raptor Changelog

## 0.10.0 -> 0.11.0
- Added the BAM and PacBio Dataset (.xml) parsing support via Pbbam. Pbbam is a big dependency, so this is turned off by default. To compile with Pbbam, type `make release-pb`.
- Refactored the sequence parsing.
- Fixed the return value when command line arguments are wrong (used to be 0, now it's 1).
- Fixed the return value when a fatal error occurs (used to be 0, no it's 1).
- Updated the build system.
- Added more unit and cram tests.
- Reorganized the code a little bit to make it more logically placed.
- Minor bugfixes.

## 0.9.14 -> 0.10.0
- Modified the composite command line options for overlapping.
- Added a `--diff` option to score mappings.
- Implemented a TravisCI build script.
