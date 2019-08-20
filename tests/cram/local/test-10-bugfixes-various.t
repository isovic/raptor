Removing the identity filter exposed the problem with calculating the Alignment Score in case of aligner_edlib. The AS was simply the number of matches, which turned out to be wrong. If only number of matches is used, then indels are ignored, and lower quality alignments can prevail
  $ ${BIN_DIR}/raptor -r ${PROJECT_DIR}/test-data/bugfixes/min_idt_issue/ctg.s1.000000F-95000-115000.fasta -g ${PROJECT_DIR}/test-data/bugfixes/min_idt_issue/all_contig.gfa2 -q ${PROJECT_DIR}/test-data/bugfixes/min_idt_issue/problematic_read.fasta -v 0 --align --out-fmt sam --bestn 0 --bestn-threshold 1.0 --min-idt 0.0
  @HD	VN:1.5
  @SQ	SN:ctg.s1.000000F:95000-115000	LN:20001
  m54081_181221_163846/11666180/49398_49554	0	ctg.s1.000000F:95000-115000	13308	46	1=1I70=1X54=2X1I1=1I1=1I1=3X18S	*	0	0	TTCATCTGCGCGGGAATGACGATTCAGAAGTTACACGAAACTCAAAAAAAACGAAACCGAACGAACCGGATTTCCGCTTTTACGGGAATGACGGCGCATAAGTTCCCGTGCGGACAGACCTAGATTCGAGAACACAAAAAAAACCAAAAAGGGGGT	*	NM:i:10	AS:i:216	QS:i:0	QE:i:138	QL:i:156	TS:i:13307	TE:i:13441	TL:i:20001	pi:i:0	pj:i:0	pn:i:1	ps:i:0
  m54081_181221_163846/11666180/49398_49554	256	ctg.s1.000000F:95000-115000	4343	46	2S2=1X5=1X21=1X1=1I1=1D5=1I1=1X8=1D12=1X15=2X45=5I1=1I2=20S	*	0	0	TTCATCTGCGCGGGAATGACGATTCAGAAGTTACACGAAACTCAAAAAAAACGAAACCGAACGAACCGGATTTCCGCTTTTACGGGAATGACGGCGCATAAGTTCCCGTGCGGACAGACCTAGATTCGAGAACACAAAAAAAACCAAAAAGGGGGT	*	NM:i:17	AS:i:178	QS:i:2	QE:i:136	QL:i:156	TS:i:4342	TE:i:4470	TL:i:20001	pi:i:2	pj:i:0	pn:i:1	ps:i:0
  m54081_181221_163846/11666180/49398_49554	256	ctg.s1.000000F:95000-115000	5151	46	1=1X2=1X27=1X1=1X6=1I1=1X6=1I1=1D12=1X15=2X14=1X25=1X4=2X1I1=1X1=2X1=1X1=2X1=2X13S	*	0	0	TTCATCTGCGCGGGAATGACGATTCAGAAGTTACACGAAACTCAAAAAAAACGAAACCGAACGAACCGGATTTCCGCTTTTACGGGAATGACGGCGCATAAGTTCCCGTGCGGACAGACCTAGATTCGAGAACACAAAAAAAACCAAAAAGGGGGT	*	NM:i:24	AS:i:144	QS:i:0	QE:i:143	QL:i:156	TS:i:5150	TE:i:5291	TL:i:20001	pi:i:3	pj:i:0	pn:i:1	ps:i:0

  $ ${BIN_DIR}/raptor -r ${PROJECT_DIR}/test-data/bugfixes/min_idt_issue/ctg.s1.000000F-95000-115000.fasta -g ${PROJECT_DIR}/test-data/bugfixes/min_idt_issue/all_contig.gfa2 -q ${PROJECT_DIR}/test-data/bugfixes/min_idt_issue/problematic_read.fasta -v 0 --align --out-fmt sam
  @HD	VN:1.5
  @SQ	SN:ctg.s1.000000F:95000-115000	LN:20001
  m54081_181221_163846/11666180/49398_49554	0	ctg.s1.000000F:95000-115000	13308	47	1=1I70=1X54=2X1I1=1I1=1I1=3X18S	*	0	0	TTCATCTGCGCGGGAATGACGATTCAGAAGTTACACGAAACTCAAAAAAAACGAAACCGAACGAACCGGATTTCCGCTTTTACGGGAATGACGGCGCATAAGTTCCCGTGCGGACAGACCTAGATTCGAGAACACAAAAAAAACCAAAAAGGGGGT	*	NM:i:10	AS:i:216	QS:i:0	QE:i:138	QL:i:156	TS:i:13307	TE:i:13441	TL:i:20001	pi:i:0	pj:i:0	pn:i:1	ps:i:0
