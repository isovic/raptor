Graph alignment on a circular reference, simplest exact mapping test case. Alignment using Edlib.
  $ ${BIN_DIR}/raptor --out-fmt paf -r ${PROJECT_DIR}/test-data/demos/demo-1-circular-synth/synth-circular.fasta.gz -g ${PROJECT_DIR}/test-data/demos/demo-1-circular-synth/graph.gfa2 -q ${PROJECT_DIR}/test-data/demos/demo-1-circular-synth/read-1-exact-match-circular.fasta.gz -v 0 --align --aligner edlib
  gi|545778205|gb|U00096.3|:40000-60000	20000	0	10000	+	synth-circ-1	100000	90000	100000	10000	10000	60	cm:i:10000	NM:i:0	AS:i:10000	pi:i:0	pj:i:0	pn:i:2	ps:i:1	cg:Z:10000=
  gi|545778205|gb|U00096.3|:40000-60000	20000	10000	20000	+	synth-circ-1	100000	0	10000	10000	10000	60	cm:i:10000	NM:i:0	AS:i:10000	pi:i:0	pj:i:1	pn:i:2	ps:i:1	cg:Z:10000=

Graph alignment on a circular reference, simplest exact mapping test case. Alignment using ksw2-single.
  $ ${BIN_DIR}/raptor --out-fmt paf -r ${PROJECT_DIR}/test-data/demos/demo-1-circular-synth/synth-circular.fasta.gz -g ${PROJECT_DIR}/test-data/demos/demo-1-circular-synth/graph.gfa2 -q ${PROJECT_DIR}/test-data/demos/demo-1-circular-synth/read-1-exact-match-circular.fasta.gz -v 0 --align --aligner ksw2-single
  gi|545778205|gb|U00096.3|:40000-60000	20000	0	10000	+	synth-circ-1	100000	90000	100000	10000	10000	60	cm:i:10000	NM:i:0	AS:i:10000	pi:i:0	pj:i:0	pn:i:2	ps:i:1	cg:Z:10000=
  gi|545778205|gb|U00096.3|:40000-60000	20000	10000	20000	+	synth-circ-1	100000	0	10000	10000	10000	60	cm:i:10000	NM:i:0	AS:i:10000	pi:i:0	pj:i:1	pn:i:2	ps:i:1	cg:Z:10000=

Graph alignment on a circular reference, simplest exact mapping test case. Alignment using ksw2-double
  $ ${BIN_DIR}/raptor --out-fmt paf -r ${PROJECT_DIR}/test-data/demos/demo-1-circular-synth/synth-circular.fasta.gz -g ${PROJECT_DIR}/test-data/demos/demo-1-circular-synth/graph.gfa2 -q ${PROJECT_DIR}/test-data/demos/demo-1-circular-synth/read-1-exact-match-circular.fasta.gz -v 0 --align --aligner ksw2-double
  gi|545778205|gb|U00096.3|:40000-60000	20000	0	10000	+	synth-circ-1	100000	90000	100000	10000	10000	60	cm:i:10000	NM:i:0	AS:i:10000	pi:i:0	pj:i:0	pn:i:2	ps:i:1	cg:Z:10000=
  gi|545778205|gb|U00096.3|:40000-60000	20000	10000	20000	+	synth-circ-1	100000	0	10000	10000	10000	60	cm:i:10000	NM:i:0	AS:i:10000	pi:i:0	pj:i:1	pn:i:2	ps:i:1	cg:Z:10000=

Another synthetic circular test case.
  $ ${BIN_DIR}/raptor --out-fmt m4 -r ${PROJECT_DIR}/test-data/graph-mapping/synth-1-rotation/ecoli-0-100000.fasta -g ${PROJECT_DIR}/test-data/graph-mapping/synth-1-rotation/graph.gfa2 -q ${PROJECT_DIR}/test-data/graph-mapping/synth-1-rotation/read-1-12kbp_rotated-fwd.fasta -v 0 --align
  read-1-12kbp_rotated-fwd gi|545778205|gb|U00096.3| -88000 100 0 0 88000 100000 0 12000 100000 100000
  read-1-12kbp_rotated-fwd gi|545778205|gb|U00096.3| -12000 100 0 88000 100000 100000 0 0 12000 100000

Yet another synthetic circular test case.
  $ ${BIN_DIR}/raptor --out-fmt m4 -r ${PROJECT_DIR}/test-data/graph-mapping/synth-1-rotation/ecoli-0-100000.fasta -g ${PROJECT_DIR}/test-data/graph-mapping/synth-1-rotation/graph.gfa2 -q ${PROJECT_DIR}/test-data/graph-mapping/synth-1-rotation/read-2-12kbp_rotated-rev.fasta -v 0 --align
  read-2-12kbp_rotated-ref gi|545778205|gb|U00096.3| -12000 100 0 0 12000 100000 1 0 12000 100000
  read-2-12kbp_rotated-ref gi|545778205|gb|U00096.3| -88000 100 0 12000 100000 100000 1 12000 100000 100000

Graph alignment of the reference sequences on a small plasmid contig set. Here, contigs are used as the target, and the contig graph as the input graph. Reference sequences are used as queries.
  $ ${BIN_DIR}/raptor --out-fmt m4 -r ${PROJECT_DIR}/test-data/graph-mapping/real-1-plasmid/asm.plasmids.fa -g ${PROJECT_DIR}/test-data/graph-mapping/real-1-plasmid/contig.s2.gfa2 -q ${PROJECT_DIR}/test-data/graph-mapping/real-1-plasmid/ref.plasmids.fa -v 0 --align
  gi|386611788|ref|NC_017637.1| ctg.s2.000000F -12567 99.7614 0 0 12593 102536 1 0 12571 102273
  gi|386611788|ref|NC_017637.1| ctg.s2.000000F -89649 99.6232 0 12593 102536 102536 1 12571 102273 102273
  gi|386607294|ref|NC_017636.1| 1 -3957 99.2172 0 0 3985 5360 0 1372 5332 5332
  gi|386607294|ref|NC_017636.1| 1 -1370 99.5627 0 3985 5360 5360 0 0 1372 5332

