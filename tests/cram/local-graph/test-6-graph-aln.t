This test uses the approximate aligner along the mapped path. It cannot pick up the small variants, because there was nothing mapped in those branches, and the distance for both branches of a SNP looks identical when mapping to the graph.
  $ ${BIN_DIR}/raptor -r ${PROJECT_DIR}/test-data/graph-aln/variants.gfa2 -g ${PROJECT_DIR}/test-data/graph-aln/variants.gfa2 -q ${PROJECT_DIR}/test-data/graph-aln/read1.fasta --align -v 0
  query1	690	0	464	+	ref1	695	0	464	464	464	60	cm:i:461	NM:i:3	AS:i:461	pi:i:0	pj:i:0	pn:i:4	ps:i:1	cg:Z:124=1X86=1X79=1X172=
  query1	690	464	469	+	sv1	5	0	5	5	5	60	cm:i:5	NM:i:0	AS:i:5	pi:i:0	pj:i:1	pn:i:4	ps:i:1	cg:Z:5=
  query1	690	469	600	+	ref1	695	464	595	131	131	60	cm:i:131	NM:i:0	AS:i:131	pi:i:0	pj:i:2	pn:i:4	ps:i:1	cg:Z:131=
  query1	690	600	690	+	ref1	695	605	695	90	90	60	cm:i:90	NM:i:0	AS:i:90	pi:i:0	pj:i:3	pn:i:4	ps:i:1	cg:Z:90=
