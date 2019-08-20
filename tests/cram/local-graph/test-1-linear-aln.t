Linear mapping on a circular reference. The following read should not be split aligned.
  $ ${BIN_DIR}/raptor --out-fmt paf -r ${PROJECT_DIR}/test-data/demos/demo-1-circular-synth/synth-circular.fasta.gz -g ${PROJECT_DIR}/test-data/demos/demo-1-circular-synth/graph.gfa2 -q ${PROJECT_DIR}/test-data/demos/demo-1-circular-synth/read-2-exact-match-linear.fasta.gz -v 0 --align --aligner edlib
  gi|545778205|gb|U00096.3|:60000-80000	20000	0	20000	+	synth-circ-1	100000	10000	30000	20000	20000	60	cm:i:20000	NM:i:0	AS:i:40000	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:20000=
