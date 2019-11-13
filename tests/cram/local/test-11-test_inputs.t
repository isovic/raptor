Run Raptor with a single RaptorDB input query file. This should work.
  $ ${BIN_DIR}/raptor-reshape -i ${PROJECT_DIR}/test-data/ecoli-small/reads.6x.fwd.fasta -o reads --block-size 0.2 --split-blocks -v 0
  > ${BIN_DIR}/raptor -r ${PROJECT_DIR}/test-data/ecoli-small/ecoli-0-100000.fasta -q reads.rdb -v 0 | wc -l | awk '{ $1=$1; print }'
  62

Run Raptor with multiple RaptorDB input query files. Only one is allowed.
  $ ${BIN_DIR}/raptor-reshape -i ${PROJECT_DIR}/test-data/ecoli-small/reads.6x.fwd.fasta -o reads --block-size 0.2 --split-blocks -v 0
  > ${BIN_DIR}/raptor -r ${PROJECT_DIR}/test-data/ecoli-small/ecoli-0-100000.fasta -q reads.rdb -q reads.rdb -v 0 2>&1 | head -n 1
  Only one RaptorDB sequence file can be loaded.

Run Raptor with multiple XML dataset input query files. Only one is allowed.
  $ ${BIN_DIR}/raptor -r ${PROJECT_DIR}/test-data/ecoli-small/ecoli-0-100000.fasta -q ${PROJECT_DIR}/test-data/small-xml/subreadset.xml -q ${PROJECT_DIR}/test-data/small-xml/subreadset.xml -v 0 2>&1 | head -n 1
  Only one XML query sequence file can be loaded.

Run Raptor with a combination of a RaptorDB and an XML input query files. This is not allowed.
  $ ${BIN_DIR}/raptor-reshape -i ${PROJECT_DIR}/test-data/ecoli-small/reads.6x.fwd.fasta -o reads --block-size 0.2 --split-blocks -v 0
  > ${BIN_DIR}/raptor -r ${PROJECT_DIR}/test-data/ecoli-small/ecoli-0-100000.fasta -q reads.rdb -q ${PROJECT_DIR}/test-data/small-xml/subreadset.xml -v 0 2>&1 | head -n 1
  A RaptorDB input cannot be specified along side to another input file.

Run Raptor with a combination of an XML and a FASTA input query files. This is not allowed.
  $ ${BIN_DIR}/raptor -r ${PROJECT_DIR}/test-data/ecoli-small/ecoli-0-100000.fasta -q ${PROJECT_DIR}/test-data/small-xml/subreadset.xml -q ${PROJECT_DIR}/test-data/ecoli-small/reads.6x.fwd.fasta -v 0 2>&1 | head -n 1
  An XML input cannot be specified along side to another input file.

Run Raptor with a combination of a RaptorDB and a FASTA input query files. This is not allowed.
  $ ${BIN_DIR}/raptor-reshape -i ${PROJECT_DIR}/test-data/ecoli-small/reads.6x.fwd.fasta -o reads --block-size 0.2 --split-blocks -v 0
  > ${BIN_DIR}/raptor -r ${PROJECT_DIR}/test-data/ecoli-small/ecoli-0-100000.fasta -q reads.rdb -q ${PROJECT_DIR}/test-data/ecoli-small/reads.6x.fwd.fasta -v 0 2>&1 | head -n 1
  A RaptorDB input cannot be specified along side to another input file.
