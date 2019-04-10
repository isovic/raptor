Graph alignment on a circular reference, simplest exact mapping test case. Alignment using Edlib.
  $ ${BIN_DIR}/raptor --out-fmt paf -r ${PROJECT_DIR}/test-data/demos/demo-1-circular-synth/synth-circular.fasta.gz -g ${PROJECT_DIR}/test-data/demos/demo-1-circular-synth/graph.gfa2 -d ${PROJECT_DIR}/test-data/demos/demo-1-circular-synth/read-1-exact-match-circular.fasta.gz -v 0 --align --aligner edlib
  gi|545778205|gb|U00096.3|:40000-60000	20000	0	10000	+	synth-circ-1	100000	90000	100000	10000	10000	60	cm:i:10000	nm:i:0	as:i:10000	pi:i:0	pj:i:0	pn:i:2	ps:i:1	cg:Z:10000=
  gi|545778205|gb|U00096.3|:40000-60000	20000	10000	20000	+	synth-circ-1	100000	0	10000	10000	10000	60	cm:i:10000	nm:i:0	as:i:10000	pi:i:0	pj:i:1	pn:i:2	ps:i:1	cg:Z:10000=

Graph alignment on a circular reference, simplest exact mapping test case. Alignment using ksw2-single.
  $ ${BIN_DIR}/raptor --out-fmt paf -r ${PROJECT_DIR}/test-data/demos/demo-1-circular-synth/synth-circular.fasta.gz -g ${PROJECT_DIR}/test-data/demos/demo-1-circular-synth/graph.gfa2 -d ${PROJECT_DIR}/test-data/demos/demo-1-circular-synth/read-1-exact-match-circular.fasta.gz -v 0 --align --aligner ksw2-single
  gi|545778205|gb|U00096.3|:40000-60000	20000	0	10000	+	synth-circ-1	100000	90000	100000	10000	10000	60	cm:i:10000	nm:i:0	as:i:10000	pi:i:0	pj:i:0	pn:i:2	ps:i:1	cg:Z:10000=
  gi|545778205|gb|U00096.3|:40000-60000	20000	10000	20000	+	synth-circ-1	100000	0	10000	10000	10000	60	cm:i:10000	nm:i:0	as:i:10000	pi:i:0	pj:i:1	pn:i:2	ps:i:1	cg:Z:10000=

Graph alignment on a circular reference, simplest exact mapping test case. Alignment using ksw2-double
  $ ${BIN_DIR}/raptor --out-fmt paf -r ${PROJECT_DIR}/test-data/demos/demo-1-circular-synth/synth-circular.fasta.gz -g ${PROJECT_DIR}/test-data/demos/demo-1-circular-synth/graph.gfa2 -d ${PROJECT_DIR}/test-data/demos/demo-1-circular-synth/read-1-exact-match-circular.fasta.gz -v 0 --align --aligner ksw2-double
  gi|545778205|gb|U00096.3|:40000-60000	20000	0	10000	+	synth-circ-1	100000	90000	100000	10000	10000	60	cm:i:10000	nm:i:0	as:i:10000	pi:i:0	pj:i:0	pn:i:2	ps:i:1	cg:Z:10000=
  gi|545778205|gb|U00096.3|:40000-60000	20000	10000	20000	+	synth-circ-1	100000	0	10000	10000	10000	60	cm:i:10000	nm:i:0	as:i:10000	pi:i:0	pj:i:1	pn:i:2	ps:i:1	cg:Z:10000=
