I ran into a test case where the same pair of query/target had 2 identical mappings reported. This should not happen.
  $ ${BIN_DIR}/raptor --out-fmt paf -r ${PROJECT_DIR}/test-data/test-no-duplicates/target.fasta -d ${PROJECT_DIR}/test-data/test-no-duplicates/query.fasta -v 0
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/1421/13892_27944	14052	1870	8139	-	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/126918/6482_30476	23994	7476	13793	398	6317	48	cm:i:35	nm:i:-1	as:i:398	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*

The same as above, but with the "overlap" mode. By design, it should no longer output multiple mappings for the same query/target pair, to make overlapping safe. Previously, it should have output the secondary mapping too. I left the original expected output here for documentation, commented out.
  $ ${BIN_DIR}/raptor --out-fmt paf -r ${PROJECT_DIR}/test-data/test-no-duplicates/target.fasta -d ${PROJECT_DIR}/test-data/test-no-duplicates/query.fasta -x ovl-raw -k 15 -w 5 -v 0
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/1421/13892_27944	14052	1870	8139	-	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/126918/6482_30476	23994	7476	13793	398	6317	48	cm:i:35	nm:i:-1	as:i:398	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
#  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/1421/13892_27944	14052	11918	13577	-	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/126918/6482_30476	23994	169	1895	140	1726	48	cm:i:13	nm:i:-1	as:i:140	pi:i:1	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
