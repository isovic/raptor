Mapping a long read to an assembly graph composed of two collapsed regions (at the front and the back) and a bubble in the middle. The bubble has 2 identical branches. The long read should span the bubble.
  $ ${BIN_DIR}/raptor -r ${PROJECT_DIR}/test-data/demos/demo-5-assembly-graph-short-bubble/ref.fasta -g ${PROJECT_DIR}/test-data/demos/demo-5-assembly-graph-short-bubble/graph.gfa2 -q ${PROJECT_DIR}/test-data/demos/demo-5-assembly-graph-short-bubble/read-1.fasta -v 0 --bestn 1
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/806/0_37037	37037	215	10833	+	part1	12000	1986	11986	5668	10000	48	cm:i:482	NM:i:-1	AS:i:5657	fg:i:0	pi:i:0	pj:i:0	pn:i:3	ps:i:1	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/806/0_37037	37037	10856	21437	+	part2b	10000	6	9939	5681	9933	48	cm:i:482	NM:i:-1	AS:i:5674	fg:i:2048	pi:i:0	pj:i:1	pn:i:3	ps:i:1	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/806/0_37037	37037	21511	36973	+	part3	27467	3	14453	7380	14450	48	cm:i:633	NM:i:-1	AS:i:7370	fg:i:2048	pi:i:0	pj:i:2	pn:i:3	ps:i:1	cg:Z:*
# This one is filtered out because only "--bestn 1" is used, so only the primary path should be printed.
# m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/806/0_37037	37037	10856	21437	+	part2a	10000	6	9939	5681	9933	48	cm:i:482	NM:i:-1	AS:i:5674	fg:i:256	pi:i:1	pj:i:0	pn:i:1	ps:i:0	cg:Z:*

  $ ${BIN_DIR}/raptor -r ${PROJECT_DIR}/test-data/demos/demo-5-assembly-graph-short-bubble/ref.fasta -g ${PROJECT_DIR}/test-data/demos/demo-5-assembly-graph-short-bubble/graph.gfa2 -q ${PROJECT_DIR}/test-data/demos/demo-5-assembly-graph-short-bubble/read-1.fasta -v 0 --bestn 10
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/806/0_37037	37037	215	10833	+	part1	12000	1986	11986	5668	10000	48	cm:i:482	NM:i:-1	AS:i:5657	fg:i:0	pi:i:0	pj:i:0	pn:i:3	ps:i:1	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/806/0_37037	37037	10856	21437	+	part2b	10000	6	9939	5681	9933	48	cm:i:482	NM:i:-1	AS:i:5674	fg:i:2048	pi:i:0	pj:i:1	pn:i:3	ps:i:1	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/806/0_37037	37037	21511	36973	+	part3	27467	3	14453	7380	14450	48	cm:i:633	NM:i:-1	AS:i:7370	fg:i:2048	pi:i:0	pj:i:2	pn:i:3	ps:i:1	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/806/0_37037	37037	10856	21437	+	part2a	10000	6	9939	5681	9933	48	cm:i:482	NM:i:-1	AS:i:5674	fg:i:256	pi:i:1	pj:i:0	pn:i:1	ps:i:0	cg:Z:*

Same as the previous test, but do not use graph-based mapping.
  $ ${BIN_DIR}/raptor -r ${PROJECT_DIR}/test-data/demos/demo-5-assembly-graph-short-bubble/ref.fasta -g ${PROJECT_DIR}/test-data/demos/demo-5-assembly-graph-short-bubble/graph.gfa2 -q ${PROJECT_DIR}/test-data/demos/demo-5-assembly-graph-short-bubble/read-1.fasta -v 0 --no-gm
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/806/0_37037	37037	215	10833	+	part1	12000	1986	11986	5668	10000	1	cm:i:482	NM:i:-1	AS:i:-1	fg:i:0	pi:i:-1	pj:i:-1	pn:i:-1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/806/0_37037	37037	10856	21437	+	part2a	10000	6	9939	5681	9933	1	cm:i:482	NM:i:-1	AS:i:-1	fg:i:0	pi:i:-1	pj:i:-1	pn:i:-1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/806/0_37037	37037	10856	21437	+	part2b	10000	6	9939	5681	9933	1	cm:i:482	NM:i:-1	AS:i:-1	fg:i:0	pi:i:-1	pj:i:-1	pn:i:-1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/806/0_37037	37037	21511	36973	+	part3	27467	3	14453	7380	14450	1	cm:i:633	NM:i:-1	AS:i:-1	fg:i:0	pi:i:-1	pj:i:-1	pn:i:-1	ps:i:0	cg:Z:*


Aligning a long read to an assembly graph composed of two collapsed regions (at the front and the back) and a bubble in the middle. The bubble has 2 identical branches. The long read should span the bubble.
  $ ${BIN_DIR}/raptor -r ${PROJECT_DIR}/test-data/demos/demo-5-assembly-graph-short-bubble/ref.fasta -g ${PROJECT_DIR}/test-data/demos/demo-5-assembly-graph-short-bubble/graph.gfa2 -q ${PROJECT_DIR}/test-data/demos/demo-5-assembly-graph-short-bubble/read-1.fasta -v 0 --align --out-fmt m4
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/806/0_37037 part1 -15914 90.2137 0 214 10847 37037 0 1986 12000 12000
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/806/0_37037 part2b -15912 90.12 0 10847 21507 37037 0 0 10000 10000
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/806/0_37037 part3 -22472 89.1247 0 21507 37037 37037 0 0 14510 27467
#  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/806/0_37037 part2a -9807 90.7317 0 10847 21507 37037 0 0 10000 10000
