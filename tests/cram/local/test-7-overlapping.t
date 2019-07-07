Test overlapping of blocks with single arcs. The "-x ovl-raw" implicitly defines "--overlap-skip-self". Blocks 0 and 0 are compared.
  $ (${BIN_DIR}/raptor-reshape -v 0 -i ${PROJECT_DIR}/test-data/ecoli-small/two_identical_long_reads.fa --block-size 0.01 -o reads --symlink) &&  ${BIN_DIR}/raptor --out-fmt paf -r reads.rdb -q reads.rdb -v 0 -x ovl-raw -k 15 -w 5 --overlap-single-arc --rdb-block-ref 0 --rdb-block-query 0

Test overlapping of blocks with single arcs. The "-x ovl-raw" implicitly defines "--overlap-skip-self". Blocks 0 and 1 are compared.
  $ (${BIN_DIR}/raptor-reshape -v 0 -i ${PROJECT_DIR}/test-data/ecoli-small/two_identical_long_reads.fa --block-size 0.01 -o reads --symlink) &&  ${BIN_DIR}/raptor --out-fmt paf -r reads.rdb -q reads.rdb -v 0 -x ovl-raw -k 15 -w 5 --overlap-single-arc --rdb-block-ref 0 --rdb-block-query 1
  read-2	24292	2	24275	+	read-1	24292	2	24275	24268	24273	60	cm:i:1756	nm:i:-1	as:i:24268	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*

Test overlapping of blocks with single arcs. The "-x ovl-raw" implicitly defines "--overlap-skip-self". Blocks 1 and 0 are compared.
  $ (${BIN_DIR}/raptor-reshape -v 0 -i ${PROJECT_DIR}/test-data/ecoli-small/two_identical_long_reads.fa --block-size 0.01 -o reads --symlink) &&  ${BIN_DIR}/raptor --out-fmt paf -r reads.rdb -q reads.rdb -v 0 -x ovl-raw -k 15 -w 5 --overlap-single-arc --rdb-block-ref 1 --rdb-block-query 0

Test overlapping of blocks with single arcs. The "-x ovl-raw" implicitly defines "--overlap-skip-self". Blocks 1 and 1 are compared.
  $ (${BIN_DIR}/raptor-reshape -v 0 -i ${PROJECT_DIR}/test-data/ecoli-small/two_identical_long_reads.fa --block-size 0.01 -o reads --symlink) &&  ${BIN_DIR}/raptor --out-fmt paf -r reads.rdb -q reads.rdb -v 0 -x ovl-raw -k 15 -w 5 --overlap-single-arc --rdb-block-ref 1 --rdb-block-query 1



Test overlapping of blocks with dual arcs. The "-x ovl-raw" implicitly defines "--overlap-skip-self". Blocks 0 and 0 are compared.
  $ (${BIN_DIR}/raptor-reshape -v 0 -i ${PROJECT_DIR}/test-data/ecoli-small/two_identical_long_reads.fa --block-size 0.01 -o reads --symlink) &&  ${BIN_DIR}/raptor --out-fmt paf -r reads.rdb -q reads.rdb -v 0 -x ovl-raw -k 15 -w 5 --rdb-block-ref 0 --rdb-block-query 0

Test overlapping of blocks with dual arcs. The "-x ovl-raw" implicitly defines "--overlap-skip-self". Blocks 0 and 1 are compared.
  $ (${BIN_DIR}/raptor-reshape -v 0 -i ${PROJECT_DIR}/test-data/ecoli-small/two_identical_long_reads.fa --block-size 0.01 -o reads --symlink) &&  ${BIN_DIR}/raptor --out-fmt paf -r reads.rdb -q reads.rdb -v 0 -x ovl-raw -k 15 -w 5 --rdb-block-ref 0 --rdb-block-query 1
  read-2	24292	2	24275	+	read-1	24292	2	24275	24268	24273	60	cm:i:1756	nm:i:-1	as:i:24268	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*

Test overlapping of blocks with dual arcs. The "-x ovl-raw" implicitly defines "--overlap-skip-self". Blocks 1 and 0 are compared.
  $ (${BIN_DIR}/raptor-reshape -v 0 -i ${PROJECT_DIR}/test-data/ecoli-small/two_identical_long_reads.fa --block-size 0.01 -o reads --symlink) &&  ${BIN_DIR}/raptor --out-fmt paf -r reads.rdb -q reads.rdb -v 0 -x ovl-raw -k 15 -w 5 --rdb-block-ref 1 --rdb-block-query 0
  read-1	24292	2	24275	+	read-2	24292	2	24275	24268	24273	60	cm:i:1756	nm:i:-1	as:i:24268	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*

Test overlapping of blocks with dual arcs. The "-x ovl-raw" implicitly defines "--overlap-skip-self". Blocks 1 and 1 are compared.
  $ (${BIN_DIR}/raptor-reshape -v 0 -i ${PROJECT_DIR}/test-data/ecoli-small/two_identical_long_reads.fa --block-size 0.01 -o reads --symlink) &&  ${BIN_DIR}/raptor --out-fmt paf -r reads.rdb -q reads.rdb -v 0 -x ovl-raw -k 15 -w 5 --rdb-block-ref 1 --rdb-block-query 1
