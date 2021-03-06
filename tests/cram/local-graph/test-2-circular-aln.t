Graph alignment on a circular reference, simplest exact mapping test case. Alignment using Edlib.
  $ ${BIN_DIR}/raptor --out-fmt paf -r ${PROJECT_DIR}/test-data/demos/demo-1-circular-synth/synth-circular.fasta.gz -g ${PROJECT_DIR}/test-data/demos/demo-1-circular-synth/graph.gfa2 -q ${PROJECT_DIR}/test-data/demos/demo-1-circular-synth/read-1-exact-match-circular.fasta.gz -v 0 --align --aligner edlib | sed -E 's/AS:i:[0-9]+[[:space:]]//g'
  gi|545778205|gb|U00096.3|:40000-60000	20000	0	10000	+	synth-circ-1	100000	90000	100000	10000	10000	60	cm:i:10000	NM:i:0	fg:i:0	pi:i:0	pj:i:0	pn:i:2	ps:i:1	cg:Z:10000=
  gi|545778205|gb|U00096.3|:40000-60000	20000	10000	20000	+	synth-circ-1	100000	0	10000	10000	10000	60	cm:i:10000	NM:i:0	fg:i:2048	pi:i:0	pj:i:1	pn:i:2	ps:i:1	cg:Z:10000=

Graph alignment on a circular reference, simplest exact mapping test case. Alignment using ksw2-single.
  $ ${BIN_DIR}/raptor --out-fmt paf -r ${PROJECT_DIR}/test-data/demos/demo-1-circular-synth/synth-circular.fasta.gz -g ${PROJECT_DIR}/test-data/demos/demo-1-circular-synth/graph.gfa2 -q ${PROJECT_DIR}/test-data/demos/demo-1-circular-synth/read-1-exact-match-circular.fasta.gz -v 0 --align --aligner ksw2-single | sed -E 's/AS:i:[0-9]+[[:space:]]//g'
  gi|545778205|gb|U00096.3|:40000-60000	20000	0	10000	+	synth-circ-1	100000	90000	100000	10000	10000	60	cm:i:10000	NM:i:0	fg:i:0	pi:i:0	pj:i:0	pn:i:2	ps:i:1	cg:Z:10000=
  gi|545778205|gb|U00096.3|:40000-60000	20000	10000	20000	+	synth-circ-1	100000	0	10000	10000	10000	60	cm:i:10000	NM:i:0	fg:i:2048	pi:i:0	pj:i:1	pn:i:2	ps:i:1	cg:Z:10000=

Graph alignment on a circular reference, simplest exact mapping test case. Alignment using ksw2-double
  $ ${BIN_DIR}/raptor --out-fmt paf -r ${PROJECT_DIR}/test-data/demos/demo-1-circular-synth/synth-circular.fasta.gz -g ${PROJECT_DIR}/test-data/demos/demo-1-circular-synth/graph.gfa2 -q ${PROJECT_DIR}/test-data/demos/demo-1-circular-synth/read-1-exact-match-circular.fasta.gz -v 0 --align --aligner ksw2-double | sed -E 's/AS:i:[0-9]+[[:space:]]//g'
  gi|545778205|gb|U00096.3|:40000-60000	20000	0	10000	+	synth-circ-1	100000	90000	100000	10000	10000	60	cm:i:10000	NM:i:0	fg:i:0	pi:i:0	pj:i:0	pn:i:2	ps:i:1	cg:Z:10000=
  gi|545778205|gb|U00096.3|:40000-60000	20000	10000	20000	+	synth-circ-1	100000	0	10000	10000	10000	60	cm:i:10000	NM:i:0	fg:i:2048	pi:i:0	pj:i:1	pn:i:2	ps:i:1	cg:Z:10000=

Another synthetic circular test case.
  $ ${BIN_DIR}/raptor --out-fmt m4 -r ${PROJECT_DIR}/test-data/graph-mapping/synth-1-rotation/ecoli-0-100000.fasta -g ${PROJECT_DIR}/test-data/graph-mapping/synth-1-rotation/graph.gfa2 -q ${PROJECT_DIR}/test-data/graph-mapping/synth-1-rotation/read-1-12kbp_rotated-fwd.fasta -v 0 --align
  read-1-12kbp_rotated-fwd gi|545778205|gb|U00096.3| -176000 100 0 0 88000 100000 0 12000 100000 100000
  read-1-12kbp_rotated-fwd gi|545778205|gb|U00096.3| -24000 100 0 88000 100000 100000 0 0 12000 100000

Yet another synthetic circular test case.
  $ ${BIN_DIR}/raptor --out-fmt m4 -r ${PROJECT_DIR}/test-data/graph-mapping/synth-1-rotation/ecoli-0-100000.fasta -g ${PROJECT_DIR}/test-data/graph-mapping/synth-1-rotation/graph.gfa2 -q ${PROJECT_DIR}/test-data/graph-mapping/synth-1-rotation/read-2-12kbp_rotated-rev.fasta -v 0 --align
  read-2-12kbp_rotated-ref gi|545778205|gb|U00096.3| -24000 100 0 0 12000 100000 1 0 12000 100000
  read-2-12kbp_rotated-ref gi|545778205|gb|U00096.3| -176000 100 0 12000 100000 100000 1 12000 100000 100000

Graph alignment of the reference sequences on a small plasmid contig set. Here, contigs are used as the target, and the contig graph as the input graph. Reference sequences are used as queries.
The 4 short alignments are all secondary alignments, and are correctly marked as secondary.
  $ ${BIN_DIR}/raptor --out-fmt m4 -r ${PROJECT_DIR}/test-data/graph-mapping/real-1-plasmid/asm.plasmids.fa -g ${PROJECT_DIR}/test-data/graph-mapping/real-1-plasmid/contig.s2.gfa2 -q ${PROJECT_DIR}/test-data/graph-mapping/real-1-plasmid/ref.plasmids.fa -v 0 --align
  gi|386611788|ref|NC_017637.1| ctg.s2.000000F -25018 99.7935 0 0 12593 102536 1 0 12571 102273
  gi|386611788|ref|NC_017637.1| ctg.s2.000000F -177984 99.6731 0 12593 102536 102536 1 12571 102273 102273
  gi|386607294|ref|NC_017636.1| 1 -7798 99.2974 0 0 3985 5360 0 1372 5332 5332
  gi|386607294|ref|NC_017636.1| 1 -2716 99.6364 0 3985 5360 5360 0 0 1372 5332
