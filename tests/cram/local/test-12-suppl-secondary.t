Run a small test sample with missing adapter. One piece should be reported as a secondary and one as primary.
  $ ${BIN_DIR}/raptor -r ${PROJECT_DIR}/test-data/secondary/test-1-adapter/ref.fasta -q ${PROJECT_DIR}/test-data/secondary/test-1-adapter/reads.fasta -v 0 | sed 's/[[:space:]]pi:.*ps:i:[0-9]//g'
  m54081_181222_131903/18219648/54727_74589	19862	18	9745	-	ctg.s1.000000F:53800-63300	9501	143	9414	4914	9271	12	cm:i:415	NM:i:-1	AS:i:4914	fg:i:16	cg:Z:*
  m54081_181222_131903/18219648/54727_74589	19862	10050	19845	+	ctg.s1.000000F:53800-63300	9501	82	9416	4841	9334	12	cm:i:408	NM:i:-1	AS:i:4841	fg:i:256	cg:Z:*

The same test example as above, but testing filtering on mapped length. The primary piece has mapped length 9727 and secondary is larger with a span of 9795. If the primary piece is filtered out by length, the secondary shouldn't be reported either.
  $ ${BIN_DIR}/raptor -r ${PROJECT_DIR}/test-data/secondary/test-1-adapter/ref.fasta -q ${PROJECT_DIR}/test-data/secondary/test-1-adapter/reads.fasta -v 0 --min-map-len 9278 | sed 's/[[:space:]]pi:.*ps:i:[0-9]//g'
  m54081_181222_131903/18219648/54727_74589	19862	18	9745	-	ctg.s1.000000F:53800-63300	9501	143	9414	4914	9271	12	cm:i:415	NM:i:-1	AS:i:4914	fg:i:16	cg:Z:*
  m54081_181222_131903/18219648/54727_74589	19862	10050	19845	+	ctg.s1.000000F:53800-63300	9501	82	9416	4841	9334	12	cm:i:408	NM:i:-1	AS:i:4841	fg:i:256	cg:Z:*

Test the bestn filter. Setting "--bestn 1" should report only the primary alignment.
  $ ${BIN_DIR}/raptor -r ${PROJECT_DIR}/test-data/secondary/test-1-adapter/ref.fasta -q ${PROJECT_DIR}/test-data/secondary/test-1-adapter/reads.fasta -v 0 --bestn 1 --bestn-threshold 1.0 | sed 's/[[:space:]]pi:.*ps:i:[0-9]//g'
  m54081_181222_131903/18219648/54727_74589	19862	18	9745	-	ctg.s1.000000F:53800-63300	9501	143	9414	4914	9271	12	cm:i:415	NM:i:-1	AS:i:4914	fg:i:16	cg:Z:*

Test the bestn filter. Setting "--bestn 10" should report at most 9 secondary alignments and 1 primary alignment. In this case, there is 1 secondary alignment.
  $ ${BIN_DIR}/raptor -r ${PROJECT_DIR}/test-data/secondary/test-1-adapter/ref.fasta -q ${PROJECT_DIR}/test-data/secondary/test-1-adapter/reads.fasta -v 0 --bestn 10 --bestn-threshold 1.0 | sed 's/[[:space:]]pi:.*ps:i:[0-9]//g'
  m54081_181222_131903/18219648/54727_74589	19862	18	9745	-	ctg.s1.000000F:53800-63300	9501	143	9414	4914	9271	12	cm:i:415	NM:i:-1	AS:i:4914	fg:i:16	cg:Z:*
  m54081_181222_131903/18219648/54727_74589	19862	10050	19845	+	ctg.s1.000000F:53800-63300	9501	82	9416	4841	9334	12	cm:i:408	NM:i:-1	AS:i:4841	fg:i:256	cg:Z:*

Test the bestn filter. Setting "--bestn 0" should report all secondary alignments and 1 primary alignment. In this case, there is 1 secondary alignment.
  $ ${BIN_DIR}/raptor -r ${PROJECT_DIR}/test-data/secondary/test-1-adapter/ref.fasta -q ${PROJECT_DIR}/test-data/secondary/test-1-adapter/reads.fasta -v 0 --bestn 0 --bestn-threshold 1.0 | sed 's/[[:space:]]pi:.*ps:i:[0-9]//g'
  m54081_181222_131903/18219648/54727_74589	19862	18	9745	-	ctg.s1.000000F:53800-63300	9501	143	9414	4914	9271	12	cm:i:415	NM:i:-1	AS:i:4914	fg:i:16	cg:Z:*
  m54081_181222_131903/18219648/54727_74589	19862	10050	19845	+	ctg.s1.000000F:53800-63300	9501	82	9416	4841	9334	12	cm:i:408	NM:i:-1	AS:i:4841	fg:i:256	cg:Z:*
