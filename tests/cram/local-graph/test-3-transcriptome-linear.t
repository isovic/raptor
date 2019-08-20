Align a synthetic spliced read to a synthetic genome + transcriptome. The synthetic read is actually spliced from real PacBio reads, and the reference is a chunk of E. Coli. The transcriptome graph annotation is fake.
  $ ${BIN_DIR}/raptor --out-fmt paf -r ${PROJECT_DIR}/test-data/demos/demo-2-transcriptome-synth/ecoli-0-100000.fasta -g ${PROJECT_DIR}/test-data/demos/demo-2-transcriptome-synth/graph.gfa2 -q ${PROJECT_DIR}/test-data/demos/demo-2-transcriptome-synth/read-1-spliced-linear.fasta -v 0
  spliced_read_1	2525	33	486	+	gi|545778205|gb|U00096.3|	100000	22060	22488	185	428	60	cm:i:18	NM:i:-1	AS:i:185	pi:i:0	pj:i:0	pn:i:4	ps:i:1	cg:Z:*
  spliced_read_1	2525	504	1183	+	gi|545778205|gb|U00096.3|	100000	40497	41133	355	636	60	cm:i:31	NM:i:-1	AS:i:355	pi:i:0	pj:i:1	pn:i:4	ps:i:1	cg:Z:*
  spliced_read_1	2525	1246	1436	+	gi|545778205|gb|U00096.3|	100000	59933	60123	95	190	60	cm:i:10	NM:i:-1	AS:i:95	pi:i:0	pj:i:2	pn:i:4	ps:i:1	cg:Z:*
  spliced_read_1	2525	1593	2145	+	gi|545778205|gb|U00096.3|	100000	83048	83579	218	531	60	cm:i:18	NM:i:-1	AS:i:218	pi:i:0	pj:i:3	pn:i:4	ps:i:1	cg:Z:*

Align a synthetic circRNA read on the transcriptome graph with a circular edge.
  $ ${BIN_DIR}/raptor -r ${PROJECT_DIR}/test-data/demos/demo-2-transcriptome-synth/ecoli-0-100000.fasta -g ${PROJECT_DIR}/test-data/demos/demo-2-transcriptome-synth/graph.gfa2 -q ${PROJECT_DIR}/test-data/demos/demo-2-transcriptome-synth/read-2-spliced-circular.fasta -v 0 --align --out-fmt m4
  spliced_read_2 gi|545778205|gb|U00096.3| -390 83.2237 0 0 323 2525 0 59887 60191 100000
  spliced_read_2 gi|545778205|gb|U00096.3| -694 65.9601 0 324 1321 2525 0 82988 83790 100000
  spliced_read_2 gi|545778205|gb|U00096.3| -704 88.3966 0 1322 1825 2525 0 22028 22502 100000
  spliced_read_2 gi|545778205|gb|U00096.3| -942 85.3881 0 1825 2525 2525 0 40493 41150 100000
