# Raptor Changelog

## 0.11.0 -> 0.12.0
- Fixed the GFA-2 output. It was missing the `$` symbol for end coordinates which match the length, and all sequences tended to be printed out for each query.
- Fixed segfault when input reference set was empty, and output format is a SAM or a BAM.
- Bugfix in src/raptor/graph_mapper.cc, where the predecessor node search was stopped too early. The code uses a binary search and a linear pass, but the break condition was wrong, stopping the loop when the end coordinate was too far away. This is resolved for now, but the best option would be to go back to using the interval tree.
- Minor updates to the command line help.
- The composite overlapping options no longer set the batch size by default (it used to be set to 1000 MB).

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
