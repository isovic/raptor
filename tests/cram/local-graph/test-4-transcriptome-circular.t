Align a synthetic circRNA read on the transcriptome graph with a circular edge.
  $ ${BIN_DIR}/raptor -r ${PROJECT_DIR}/test-data/demos/demo-2-transcriptome-synth/ecoli-0-100000.fasta -g ${PROJECT_DIR}/test-data/demos/demo-2-transcriptome-synth/graph.gfa2 -q ${PROJECT_DIR}/test-data/demos/demo-2-transcriptome-synth/read-2-spliced-circular.fasta -v 0 --align --out-fmt m4
  spliced_read_2 gi|545778205|gb|U00096.3| -390 87.9257 0 0 323 2525 0 59887 60191 100000
  spliced_read_2 gi|545778205|gb|U00096.3| -694 75.4263 0 324 1321 2525 0 82988 83790 100000
  spliced_read_2 gi|545778205|gb|U00096.3| -704 91.0537 0 1322 1825 2525 0 22028 22502 100000
  spliced_read_2 gi|545778205|gb|U00096.3| -942 89.7143 0 1825 2525 2525 0 40493 41150 100000
