# raptor
Graph-based mapping of long sequences, noisy or HiFi

Master: [![Build Status](https://travis-ci.com/isovic/raptor.svg?token=RuhN6p8VX1u4RgANMqMt&branch=master)](https://travis-ci.com/isovic/raptor) Dev: [![Build Status](https://travis-ci.com/isovic/raptor.svg?token=RuhN6p8VX1u4RgANMqMt&branch=dev)](https://travis-ci.com/isovic/raptor)

## Contents
- [Features and Overview](#features)
- [Work in progress](#wip)
- [Disclaimer](#disclaimer)
- [Building](#building)
  - [Building with BAM support (both input and output)](#building_with_bam)
  - [Building without the BAM support)](#building_no_bam)
  - [Dependencies)](#dependencies)
- [Example usage](#example_usage)
- [Graph mapping](#graph_mapping)
- [Output formats](#output_formats)
- [Algorithmic overview](#algorithmic_overview)
- [License](#license)

## <a name="features"></a>Features and Overview
Raptor is a very versatile and fast **graph based sequence mapper/aligner** with a large number of features and more to come:
- Sequence-to-Graph mapping and path alignment.
- Linear sequence-to-sequence mapping/alignment when no graph is provided.
- Overlapping of sequences.
- Very fast mapping - on the order of magnitude of Minimap2, though somewhat slower.
- Minimizer indexing.
- Homopolymer suppression for mapping.
- Large number of input formats supported: `fasta`, `fastq`, `gfa1`, `gfa2`, `sam`, `bam`, `rdb` (RaptorDB), and `gzipped` versions of those. Also, PacBio `.xml`.
- A number of supported output formats: `paf`, `sam`, `bam`, `m4`, `mhap`, `gfa2`.
- Interchangeable aligners. For now: Edlib, KSW2; easy to add new ones.
- RaptorDB - a format similar to `.fai`, but supports an unlimited number of input files, in all supported formats listed above. Allows random access to any sequence or blocks of sequences (e.g. useful for overlapping on an HPC). Can be constructed via `raptor-reshape`.
- Can be used as a library. The Meson build is configured to produce a shared library automatically.
- Very large portions of code are covered with unit and Cram tests.

Note that it's potentially possible to also perform overlapping of sequences with graphs, but this hasn't yet been tested.

## <a name="wip"></a>Work in progress
Features:
- **Graph alignment**, to support alignment over small nucleotide variations. (Current graph mapping implementation is better suited for graphs with larger variations, such as assembly graphs, structural variations, transcriptomes, and similar).
- Mapping with inversions.
- RNA-seq/Isoseq mapping.
- Graph-based sequence simulator.

Thorough benchmarking is still on the to do list:
- Comparison with existing linear mappers.
- Comparison with existing graph-based mappers/aligners.
- Resource consumption.

## <a name="disclaimer"></a>Disclaimer
This is a proof-of-concept work in further rapid development of features.

Graph-based mapping in `raptor` is currently best suited for graphs representing large variations such as:
- Assembly graphs (e.g. `p_ctg` and `a_ctg` in Falcon)
- Haplotig graphs (e.g. from Falcon-Unzip)
- Reference graphs representing structural variations
- Circular genomes
- Transcriptomes
- (and similar)

At the moment, `raptor` does not yet align through small nucleotide variations (SNPs and short indels). **This feature is currently in development.**

## <a name="building"></a>Building
Raptor uses Meson as the build system, but the entire process is wrapped using Makefile to simplify the command line usage.

### <a name="building_with_bam"></a>Building with BAM support (both input and output)
This will compile `raptor` with support for BAM files, as well as the PacBio Dataset format (`.xml` files).
BAM support is useful, for example, to polish contig sequences using `arrow` (`raptor` will automatically parse and apply any PacBio specific tags/headers).
To build, simply type:  
```
make release-pb
```

To test the installation, run:
```
make unit
make cram
```
This option will require the `pbbam` library which is a relatively big dependency with more than a hundred source files. The `pbbam` in itself is already included as a subproject, but it also has several of it's own dependencies:
- Htslib - also wrapped with Meson's subprojects.
- Boost - wrapping is not posibble for this one, so it either requires manual setup (see `.travis.yml`) or a system dependency.

To compile without the `.bam`/`.xml` support, please read the next subsection.

### <a name="building_no_bam"></a>Building without the BAM support
To build without `pbbam`, type:
```
make release
```
This will still support `.sam` output.

### <a name="dependencies"></a>Dependencies
- Meson 0.48 (to install: `pip3 install meson`)
- C++14 (either GCC >= 6.4 or Clang >= 4.0)
- Optionally, if compiling with BAM support: Boost 1.67. `pbbam` and `htslib` are wrapped within the subprojects.

## <a name="example_usage"></a>Example usage

1. Basic run with the default parameters. This will generate mappings in the .paf format, without alignment.
```
raptor -r ref.fa -q reads.fasta > out.paf
raptor -r ref.fa -q reads.fasta -o out.paf
```

2. Mapping to an assembly graph:
If the graph is asymmetric (reverse complement sequences of each edge is _not_ provided), use the `-g` option. This applies to almost all assembly graphs:
```
raptor -r ctg.fa -g ctg.gfa -q reads.fasta -o out.paf
```
For symmetric graphs, use the capital `-G`, and the symmetric edges will not be added internally:
```
raptor -r ctg.fa -G ctg.gfa -q reads.fasta -o out.paf
```

3. Producing aligned `paf` output:
```
raptor --align -r ref.fa -q reads.fasta -o out.paf
```

4. Mapping to an assembly graph with the alignment turned on:
```
raptor --align  -r all_ctg.fa -g all_ctg.gfa -q reads.fasta -o out.paf
```

5. Producing aligned `sam` output:
```
raptor --align --out-fmt sam -r ref.fa -q reads.fasta -o out.sam
```

6. Map only to a specific region of a reference sequence:
```
raptor -r ref.fa -q reads.fasta --region "000000F:0-10000" -o out.paf
```

7. Hompolymer suppression enables higher sensitivity:
```
raptor -r ref.fa -q reads.fasta -o out.paf --hp-suppress
```

8. Combining HP-suppression and larger `k-mer` and window size:
```
raptor -r ref.fa -q reads.fasta -o out.paf --hp-suppress -k 19 -w 10
```

9. Other output formats are also supported in both the aligned and unalignead mode:
Without alignments:
```
raptor --out-fmt sam -r ref.fa -q reads.fasta -o out.sam
raptor --out-fmt bam -r ref.fa -q reads.fasta -o out.bam
raptor --out-fmt paf -r ref.fa -q reads.fasta -o out.paf
raptor --out-fmt mhap -r ref.fa -q reads.fasta -o out.mhap
raptor --out-fmt m4 -r ref.fa -q reads.fasta -o out.m4
raptor --out-fmt gfa2 -r ref.fa -q reads.fasta -o out.gfa2
```
With alignments:
```
raptor --align --out-fmt sam -r ref.fa -q reads.fasta -o out.sam
raptor --align --out-fmt bam -r ref.fa -q reads.fasta -o out.bam
raptor --align --out-fmt paf -r ref.fa -q reads.fasta -o out.paf
raptor --align --out-fmt m4 -r ref.fa -q reads.fasta -o out.m4
raptor --align --out-fmt gfa2 -r ref.fa -q reads.fasta -o out.gfa2
```

10. Output results to stdout:
```
raptor -k 19 -w 10 -r ref.fa -q reads.fasta
```

11. Sequences can be in several input formats, including `fasta`, `fastq`, `gfa1`, `gfa2`, `gfa` (autodetect version), `sam`, `bam`, `rdb` (RaptorDB), and `gzipped` versions of those:
```
raptor -r layout.gfa -q reads.fastq
```
The query sequences can also be in PacBio's `.xml` dataset format.

12. Controlling the amount of mappings/alignments which get output.  There are two options which control the amount of secondary/supplementary alignments: `--bestn` and `--bestn-threshold`.
If `--bestn <N>` is specified, then at most `N` highest scoring mappings/alignments will be output.
If `N <= 0`, then `--bestn-threshold <P>` is taken into account, and all mappings with the score `>= (1.0 - P) * max_score` will be output.
```
raptor -r ref.fa -q reads.fasta -o out.paf --bestn 1
raptor -r ref.fa -q reads.fasta -o out.paf --bestn 0 --bestn-threshold 0.05
```

13. Filtering the output by mapping quality and/or minimum percent identity:
```
# Works for both the mapping and alignment modes.
raptor -r ref.fa -q reads.fasta -o out.paf --mapq 3

# Functional only for alignment (needed to calculate the identity)
raptor -r ref.fa -q reads.fasta -o out.paf --align --min-idt 75.0
```

14. Prevent alignment extension in the alignment mode. Using the extensions might slow things down, but the alignments are usually longer.
```
raptor -r ref.fa -q reads.fasta -o out.paf --align --no-ext
```

15. Features useful for chunking and distributed running:
```
# Each of these could be qsubbed.
raptor -r ref.fa -q reads.fasta -o out.paf --start 0 --num-reads 100000
raptor -r ref.fa -q reads.fasta -o out.paf --start 100000 --num-reads 100000
raptor -r ref.fa -q reads.fasta -o out.paf --start 200000 --num-reads 100000
raptor -r ref.fa -q reads.fasta -o out.paf --start 300000 --num-reads 100000
raptor -r ref.fa -q reads.fasta -o out.paf --start 500000 --num-reads 100000
```

16. Specify multiple input reference files.
```
raptor -r ref1.fa -r ref2.fa -r ref3.fa -q reads.fasta -o out.paf
```

17. Specify multiple input query files.
```
raptor -r ref.fa -q reads1.fasta -q reads2.fasta -q reads3.fasta -o out.paf
```

18. Overlapping.
Overlapping low-accuracy reads:
```
raptor -x ovl-raw -r reads.fasta -q reads.fasta -o out.paf
```

Overlapping high-accuracy reads (CCS or preads):
```
raptor -x ovl-hifi -r reads.fasta -q reads.fasta -o out.paf
```

Align each overlap:
```
raptor -x ovl-raw --align -r reads.fasta -q reads.fasta -o out.paf
```

Any format is supported as well:
```
raptor -x ovl-raw --align -r reads.fasta -q reads.fasta --out-fmt bam -o out.bam
```

19. RaptorDB - Random access lookup DB.
RaptorDB is a format similar to `.fai`, but it's more generic. It can index multiple input files in any of the input formats, and enable random access to any particular sequence in the database, as well as partition the sequences into blocks.
To create a RaptorDB, a tool called `raptor-reshape` is used.
The produced `.rdb` file is a plain text file, which can easily be edited to filter the list of input sequences (e.g. filtering reads by some parameter).

Creating a RaptorDB from multiple input files using symlinking (sequences will not be copied), and split it into blocks of 400 MB:
```
raptor-reshape -i reads.fasta --out-prefix out --block-size 400 --symlink
```
Multiple inputs can be specified. Format is determined automatically based on the extension:
```
raptor-reshape -i reads1.fasta -i reads2.fastq -i reads3.bam -i reads4.gfa1 --out-prefix out --block-size 400 --symlink
```

Creating a RaptorDB from multiple input files, and copy the sequences into a new FASTA file. RaptorDB will be split into blocks of 400MB:
```
raptor-reshape -i subreads1.fasta -i subreads2.bam --out-fmt fasta --out-prefix out --block-size 400
```

Each block can also be written into a separate FASTA file for convenience (e.g. distributing jobs on a cluster):
```
raptor-reshape -i subreads1.fasta -i subreads2.bam --out-fmt fasta --out-prefix out --block-size 400 --split-blocks
```

## <a name="graph_mapping"></a>Graph Mapping

`raptor` can currently read both GFA-1 and GFA-2 formats to define the graph. The version of the GFA is determined automatically from the GFA header.

Some notes:
1. The graph should either be `symmetric` or `asymmetric`.

2. A `symmetric` graph is when all edges of the graph are explicitly defined. For example, this refers to both strands of the same edge (forward and reverse complement). In this case, defining a graph to map to a circular reference would need an edge for both the forward and the reverse strand, i.e.:
```
H	VN:Z:2.0
S	gi|545778205|gb|U00096.3|	4641652	*
E	edge1	gi|545778205|gb|U00096.3|+	gi|545778205|gb|U00096.3|+	4641652$	4641652$	0	0	*
E	edge1	gi|545778205|gb|U00096.3|-	gi|545778205|gb|U00096.3|-	4641652$	4641652$	0	0	*
```
The `symmetric` graphs can be loaded using the `-G` command line argument.

3. An `asymmetric` graph is the one in which the opposite strand of each edge is implicit, and `raptor` fills out the symmetric arcs automatically. This can greatly simplify the definition of the graph, and can enable straightforward mapping to assembly graphs. For example, the same circular reference could simply be defined using:
```
H	VN:Z:2.0
S	gi|545778205|gb|U00096.3|	4641652	*
E	edge1	gi|545778205|gb|U00096.3|+	gi|545778205|gb|U00096.3|+	4641652$	4641652$	0	0	*
```
The `asymmetric` graphs can be loaded using the `-g` command line argument.

4. The GFA-1/GFA-2 formats define coordinates to always be on the fwd strand of the sequence. Keep this in mind when designing edges which enter/exit the rev strand of a sequence.

## <a name="output_formats"></a>Output formats
Raptor supports multiple output formats, concretely: `paf`, `sam`, `bam`, `m4`, `mhap`, `gfa2`.

Each of these formats can be used with any combination of parameter options (even with or without alignment), but only some encode the graph-based mapping information. These are: `paf`, `sam`, `bam` and `gfa2`.

The graph-based mapping information is encoded in the form of SAM tags:
- `pi:i` - Path ID.
- `pj:i` - Segment in the path.
- `pn:i` - Number of segments in the path.
- `ps:i` - This is `1` if the path is split-aligned via graph edges, otherwise it's `0`.

This is the current definition, which may likely change as `raptor` evolves.
Concretely, it is likely that the graph information will be encoded into tags in the future, such as an adjacency list, to allow for full graph reconstruction.

## <a name="algorithmic_overview"></a>Algorithmic overview

1. Input is a set of target sequences, a set of query sequences and an optional graph in one of the supported formats.
2. Target sequences are indexed independently.
3. A `SegmentGraph` is constructed from the input graph file, where nodes are the target sequences.
4. For a query, all k-mer hits are looked-up and chained into anchors linearly (using the common dynamic programming approach).
5. Anchors are broken on any branching point in each target sequence.
6. An `AnchorGraph` is constructed, where nodes are anchors, and edges are either: (1) explicit edges from the SegmentGraph, or (2) implicit edges between neighboring anchors on the same target sequence. This allows initial linear anchoring to be more stringent in bandwidth.
7. Graph-based dynamic programming is applied on the `AnchorGraph` to chain the anchors over the graph.
8. Best scoring paths are extracted.
9. Paths are aligned (optionally).

## <a name="license"></a>License
Raptor is licensed under BSD 3-Clause Clear License.

Credit needs to be given to authors of libraries which Raptor utilizes for various purposes. These are:
It uses some of the open source libraries internally. These include:
- Edlib - Author: Martin Sosic, License: "MIT", Link: https://github.com/Martinsos/edlib
- KSW2 and Kseq - Author: Heng Li, License: "MIT", Link: https://github.com/lh3/minimap2/blob/master/ksw2.h
- ThreadPool - Author: Robert Vaser, License: "MIT", Link: https://github.com/rvaser/thread_pool
- IntervalTree - Author: Erik Garrison, License: "MIT", Link: https://github.com/ekg/intervaltree
- SparseHash - Author: Google Inc., "License BSD 3-Clause "New" or "Revised"", Link: https://github.com/sparsehash/sparsehash

## Author
Ivan Sovic, 2017-present.
