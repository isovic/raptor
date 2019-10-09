# Raptor Changelog

## 0.18.3 -> 0.19.0.
- Modified how the `--min-map-len` option was applied. Previously, even if alignment extension was used, `--min-map-len` was still applied before alignment (at the mapping stage), which means that even if the mapping was extended, it would still be filtered out. Now, the value of `min_map_len / 2` is used to filter the mappings in case the extension will be applied, and then the full `min_map_len` is applied after alignment. If alignment/flank extension is not applied, then the full `min_map_len` is used right away.
- Abstracted the index.
- Implemented a new index option - "dense" indexing which doesn't use minimizers. This is available via the `--index-type` command line option.

## 0.18.2 -> 0.18.3.
- Further significant reduction in memory consumption. All intermediate information (linear mapping, graph mapping and alignment) can be cleared after the query has been processed, because only the regions are important. For larger batches and more accurate data, this provides a significant reduction in required memory.
- Refactoring of the RaptorResults.

## 0.18.1 -> 0.18.2.
- Much lower memory consumption. Previously, the seed hits were never released after mapping was performed, which means that for a batch of sequences they would accumulate until the batch was done. The seed hits are no longer required after all mapping and alignment stages are complete, so the MappingWorker_ now clears the target hits vector. This is optional, in case the API user wants to retain the hits.

## 0.18.0 -> 0.18.1
- Refactored the index, and abstracted the index interface. Created a factory for the index (which currently only has a MinimizerIndex option).
- Static library linking now.

## 0.17.0 -> 0.18.0
- MinimizerIndex was revamped, and is now faster than before (~30% on a 1GB reference).
- Makefile now has a 'debug2' rule which compiles with debug symbols but witohut ASAN.

## 0.16.5 -> 0.17.0
- Fixed a bug which truncated the read groups from the header, and kept only one read group which propagated into the output. Added Cram tests to verify the solution.
- Added unit tests for the new ScoreCigarAlignment.
- Added a feature ("--no-gm") which allows to skip graph-mapping. This can add a little bit of speed for some contexts, e.g. overlapping. Hasn't yet been thoroughly tested.
- Updated the GraphSim and added a simple error-rate model.
- Added a tool to evaluate simulated results (graphsim-eval.py).
- Updated the time2 rule in the Makefile.
- Added a log output of the number of threads at the beginning of RunRaptor.

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
