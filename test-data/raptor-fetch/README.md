```
${BIN_DIR}/raptor-reshape -i ${PROJECT_DIR}/test-data/ecoli-small/reads.6x.fwd.fasta -o reads --symlink -v 0
raptor -x ovl-raw -r reads.rdb -q reads.rdb -o out.m4 --out-fmt m4
head -n 39 out.m4 > small_input.m4
```

