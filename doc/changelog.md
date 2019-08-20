# Raptor Changelog

## 0.16.4 -> 0.16.5
- Now explicitly failing if a PacBio .xml file is not a subreadset or an alignmentset.
- Added some Cram tests.
- The alignment score (AS) for the edit distance based alignment is now calculated by rescoring the alignment path (CIGAR vector). Previously, it was only the number of matches. This was not a good choice in case an alignment was riddled with indels, because the AS would still be high, and the final sorting is performed on the AS. Now, this closely reflects what a non-edit-distance-based aligners would produce.
- Updated the Cram tests according to the new AS calculation. Added a new Cram test file which will house tests for bugfixes. The test case which exposes the AS issue was added into that test file.
- Renamed the piecewise penalties from `(v, u)` to `(open, ext)` to reduce potential future confusion.

## 0.16.3 -> 0.16.4
- Fixed a segfault when a PacBio XML was used as the reference input to Raptor.
- Minor change in debug output (now writing the number of sequences in the index to screen).

## 0.16.2 -> 0.16.3
- No longer failing with a fatal error when the input RaptorDB is empty, and a block is specified.
- Raptor-reshape in the --symlink mode now first outputs the list of files, and then the sequences. This is important because some BAM files might be empty, and if the files aren't listed then the headers will not be parsed and forwarded downstream to the output BAM.

## 0.16.1 -> 0.16.2
- Updated the meson.build to support the latest Pbbam/Pbcopper updates.

## 0.16.0 -> 0.16.1
- Added a Meson option `tests` to optionally compile with or without unit tests.
- New Makefile rule `release-pb-no-tests` which compiles without unit tests.
- More Cram tests for the `raptor-fetch` tool.
- Refactoring of `raptor-fetch` and resolving warnings.

## 0.15.0 -> 0.16.0
- Added `raptor-fetch` to allow fetching of sequences from the RaptorDB. It supports 3 modes: `fetch` for plain fetching of sequences, `clip` which extracts subregions of the sequences based on the input BED file, and `erc` which creates a pileup of sequences for error correction (using fc_consensus) based on the input M4 overlap file.

## 0.14.0 -> 0.15.0
- The default batch size is now 500MB instead of 200MB.
- Updated the .travis.yml to install Samtools so that BAM Cram tests don't fail.
- Minor change in logging.
- Updated the README.

## 0.13.0 -> 0.14.0
- The `nm` tag in the SAM/PAF was now renamed to `NM` to comply with the standard.
- Also, changed `as:i:` -> `AS:i:`.
- Fixed the flag in BAM output. It didn't mark secondary/supplementary alignments.
- Added tests for SAM/BAM flag.
- Updated cram tests to match the new tags.
- Added a default `out-prefix` to `raptor-reshape`.
- Added the custom tags to the BAM output.

## 0.12.0 -> 0.13.0
- Important change in command line usage. The `-d`/`--reads` command line argument is now changed with `-q`/`--query`.

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
