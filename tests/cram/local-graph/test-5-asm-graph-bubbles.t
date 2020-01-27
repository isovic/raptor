Mapping a long read to an assembly graph composed of two collapsed regions (at the front and the back) and a bubble in the middle. The bubble has 2 identical branches. The long read should span the bubble.
  $ ${BIN_DIR}/raptor -r ${PROJECT_DIR}/test-data/demos/demo-5-assembly-graph-short-bubble/ref.fasta -g ${PROJECT_DIR}/test-data/demos/demo-5-assembly-graph-short-bubble/graph.gfa2 -q ${PROJECT_DIR}/test-data/demos/demo-5-assembly-graph-short-bubble/read-1.fasta -v 0 --bestn 1
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/806/0_37037	37037	215	10833	+	part1	12000	1986	11986	5668	10000	60	cm:i:482	NM:i:-1	AS:i:5657	fg:i:0	pi:i:0	pj:i:0	pn:i:3	ps:i:1	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/806/0_37037	37037	10856	21437	+	part2b	10000	6	9939	5681	9933	3	cm:i:482	NM:i:-1	AS:i:5674	fg:i:2048	pi:i:0	pj:i:1	pn:i:3	ps:i:1	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/806/0_37037	37037	21511	36973	+	part3	27467	3	14453	7380	14450	60	cm:i:633	NM:i:-1	AS:i:7370	fg:i:2048	pi:i:0	pj:i:2	pn:i:3	ps:i:1	cg:Z:*
# This one is filtered out because only "--bestn 1" is used, so only the primary path should be printed.
# m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/806/0_37037	37037	10856	21437	+	part2a	10000	6	9939	5681	9933	3	cm:i:482	NM:i:-1	AS:i:5674	fg:i:256	pi:i:1	pj:i:0	pn:i:1	ps:i:0	cg:Z:*

  $ ${BIN_DIR}/raptor -r ${PROJECT_DIR}/test-data/demos/demo-5-assembly-graph-short-bubble/ref.fasta -g ${PROJECT_DIR}/test-data/demos/demo-5-assembly-graph-short-bubble/graph.gfa2 -q ${PROJECT_DIR}/test-data/demos/demo-5-assembly-graph-short-bubble/read-1.fasta -v 0 --bestn 10
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/806/0_37037	37037	215	10833	+	part1	12000	1986	11986	5668	10000	60	cm:i:482	NM:i:-1	AS:i:5657	fg:i:0	pi:i:0	pj:i:0	pn:i:3	ps:i:1	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/806/0_37037	37037	10856	21437	+	part2b	10000	6	9939	5681	9933	3	cm:i:482	NM:i:-1	AS:i:5674	fg:i:2048	pi:i:0	pj:i:1	pn:i:3	ps:i:1	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/806/0_37037	37037	21511	36973	+	part3	27467	3	14453	7380	14450	60	cm:i:633	NM:i:-1	AS:i:7370	fg:i:2048	pi:i:0	pj:i:2	pn:i:3	ps:i:1	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/806/0_37037	37037	10856	21437	+	part2a	10000	6	9939	5681	9933	3	cm:i:482	NM:i:-1	AS:i:5674	fg:i:256	pi:i:1	pj:i:0	pn:i:1	ps:i:0	cg:Z:*

Same as the previous test, but do not use graph-based mapping.
  $ ${BIN_DIR}/raptor -r ${PROJECT_DIR}/test-data/demos/demo-5-assembly-graph-short-bubble/ref.fasta -g ${PROJECT_DIR}/test-data/demos/demo-5-assembly-graph-short-bubble/graph.gfa2 -q ${PROJECT_DIR}/test-data/demos/demo-5-assembly-graph-short-bubble/read-1.fasta -v 0 --no-gm
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/806/0_37037	37037	21511	36973	+	part3	27467	3	14453	7380	14450	60	cm:i:633	NM:i:-1	AS:i:-1	fg:i:0	pi:i:-1	pj:i:-1	pn:i:-1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/806/0_37037	37037	215	10833	+	part1	12000	1986	11986	5668	10000	60	cm:i:482	NM:i:-1	AS:i:-1	fg:i:2048	pi:i:-1	pj:i:-1	pn:i:-1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/806/0_37037	37037	10856	21437	+	part2b	10000	6	9939	5681	9933	3	cm:i:482	NM:i:-1	AS:i:-1	fg:i:2048	pi:i:-1	pj:i:-1	pn:i:-1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/806/0_37037	37037	10856	21437	+	part2a	10000	6	9939	5681	9933	3	cm:i:482	NM:i:-1	AS:i:-1	fg:i:256	pi:i:-1	pj:i:-1	pn:i:-1	ps:i:0	cg:Z:*


Aligning a long read to an assembly graph composed of two collapsed regions (at the front and the back) and a bubble in the middle. The bubble has 2 identical branches. The long read should span the bubble.
The second branch of the bubble is reported as a secondary alignment.
  $ ${BIN_DIR}/raptor -r ${PROJECT_DIR}/test-data/demos/demo-5-assembly-graph-short-bubble/ref.fasta -g ${PROJECT_DIR}/test-data/demos/demo-5-assembly-graph-short-bubble/graph.gfa2 -q ${PROJECT_DIR}/test-data/demos/demo-5-assembly-graph-short-bubble/read-1.fasta -v 0 --align --out-fmt paf | sed 's/[[:space:]]cg:Z:.*//g'
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/806/0_37037	37037	214	10847	+	part1	12000	1986	12000	10633	10014	60	cm:i:9797	NM:i:980	AS:i:15914	fg:i:0	pi:i:0	pj:i:0	pn:i:3	ps:i:1
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/806/0_37037	37037	10847	21507	+	part2b	10000	0	10000	10660	10000	3	cm:i:9807	NM:i:988	AS:i:15912	fg:i:2048	pi:i:0	pj:i:1	pn:i:3	ps:i:1
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/806/0_37037	37037	21507	37037	+	part3	27467	0	14510	15530	14510	60	cm:i:14182	NM:i:1578	AS:i:22472	fg:i:2048	pi:i:0	pj:i:2	pn:i:3	ps:i:1
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/806/0_37037	37037	10855	21507	+	part2a	10000	5	10000	10652	9995	3	cm:i:9802	NM:i:985	AS:i:15912	fg:i:256	pi:i:1	pj:i:0	pn:i:1	ps:i:0

Aligning a long read to an assembly graph composed of two collapsed regions (at the front and the back) and a bubble in the middle. The bubble has 2 identical branches. The long read should span the bubble.
Only the primary alignment is supposed to be reported (because of the "--bestn 1" option).
  $ ${BIN_DIR}/raptor -r ${PROJECT_DIR}/test-data/demos/demo-5-assembly-graph-short-bubble/ref.fasta -g ${PROJECT_DIR}/test-data/demos/demo-5-assembly-graph-short-bubble/graph.gfa2 -q ${PROJECT_DIR}/test-data/demos/demo-5-assembly-graph-short-bubble/read-1.fasta -v 0 --align --out-fmt paf --bestn 1 | sed 's/[[:space:]]cg:Z:.*//g'
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/806/0_37037	37037	214	10847	+	part1	12000	1986	12000	10633	10014	60	cm:i:9797	NM:i:980	AS:i:15914	fg:i:0	pi:i:0	pj:i:0	pn:i:3	ps:i:1
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/806/0_37037	37037	10847	21507	+	part2b	10000	0	10000	10660	10000	3	cm:i:9807	NM:i:988	AS:i:15912	fg:i:2048	pi:i:0	pj:i:1	pn:i:3	ps:i:1
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/806/0_37037	37037	21507	37037	+	part3	27467	0	14510	15530	14510	60	cm:i:14182	NM:i:1578	AS:i:22472	fg:i:2048	pi:i:0	pj:i:2	pn:i:3	ps:i:1
