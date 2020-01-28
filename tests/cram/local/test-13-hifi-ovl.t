Dovetail overlap 5prime test.
  $ ${BIN_DIR}/raptor -x ovl-hifi -r ${PROJECT_DIR}/test-data/hifi-ovl/test-1-dovetail/reads.pile1-5prime.fasta -q ${PROJECT_DIR}/test-data/hifi-ovl/test-1-dovetail/reads.pile1-5prime.fasta -k 30 -w 80 -v 0 -t 1 -n 1 --out-fmt m4
  m64030_190330_071939/101844710/ccs m64030_190330_071939/60686570/ccs -4972 98.9535 0 0 9150 11811 0 981 10059 10059
  m64030_190330_071939/101844710/ccs m64030_190330_071939/62523930/ccs -4592 99.9367 0 0 7895 11811 1 0 7894 11307
  m64030_190330_071939/101844710/ccs m64030_190330_071939/61737670/ccs -3342 99.8175 0 0 6031 11811 0 2564 8590 8590
  m64030_190330_071939/101844710/ccs m64030_190330_071939/60884410/ccs -2822 98.9714 0 0 5007 11811 0 4181 9139 9139
  m64030_190330_071939/101844710/ccs m64030_190330_071939/61276980/ccs -2124 99.7655 0 0 3845 11811 1 0 3838 10150

Dovetail overlap 3prime test.
  $ ${BIN_DIR}/raptor -x ovl-hifi -r ${PROJECT_DIR}/test-data/hifi-ovl/test-1-dovetail/reads.pile2-3prime.fasta -q ${PROJECT_DIR}/test-data/hifi-ovl/test-1-dovetail/reads.pile2-3prime.fasta -k 30 -w 80 -v 0 -t 1 -n 1 --out-fmt m4
  m64030_190330_071939/101909220/ccs m64030_190330_071939/165937510/ccs -4352 99.9056 0 2111 9530 9530 1 2475 9898 9898
  m64030_190330_071939/101909220/ccs m64030_190330_071939/22611030/ccs -2538 99.1535 0 4096 9530 9530 0 0 5435 9487
  m64030_190330_071939/101909220/ccs m64030_190330_071939/25626700/ccs -1902 99.9086 0 6247 9530 9530 0 0 3286 10208
  m64030_190330_071939/101909220/ccs m64030_190330_071939/154862090/ccs -1664 99.5862 0 6629 9530 9530 1 5043 7943 8044
  m64030_190330_071939/101909220/ccs m64030_190330_071939/157549180/ccs -757 99.0087 0 7809 9530 9530 1 8301 10016 10128

Dovetail overlap all fwd oriented.
  $ ${BIN_DIR}/raptor -x ovl-hifi -r ${PROJECT_DIR}/test-data/hifi-ovl/test-1-dovetail/reads.pile3-fwd.fasta -q ${PROJECT_DIR}/test-data/hifi-ovl/test-1-dovetail/reads.pile3-fwd.fasta -k 30 -w 80 -v 0 -t 1 -n 1 --out-fmt m4
  m64030_190330_071939/102303370/ccs m64030_190330_071939/106038710/ccs -4060 99.9105 0 1793 8496 8602 0 0 6707 8913
  m64030_190330_071939/102303370/ccs m64030_190330_071939/109642940/ccs -4003 99.8513 0 0 6723 8602 0 1987 8716 8823
  m64030_190330_071939/102303370/ccs m64030_190330_071939/108135170/ccs -3324 99.1565 0 0 7006 8602 0 656 7651 7651
  m64030_190330_071939/102303370/ccs m64030_190330_071939/114034050/ccs -955 99.7762 0 0 1787 8602 0 7139 8928 8928
  m64030_190330_071939/102303370/ccs m64030_190330_071939/10946880/ccs -712 99.918 0 7276 8496 8602 0 0 1221 9090

Dovetail overlap all rev oriented.
  $ ${BIN_DIR}/raptor -x ovl-hifi -r ${PROJECT_DIR}/test-data/hifi-ovl/test-1-dovetail/reads.pile4-rev.fasta -q ${PROJECT_DIR}/test-data/hifi-ovl/test-1-dovetail/reads.pile4-rev.fasta -k 30 -w 80 -v 0 -t 1 -n 1 --out-fmt m4
  m64030_190330_071939/102172020/ccs m64030_190330_071939/43124430/ccs -5719 99.7056 0 0 9869 10635 1 0 9851 10977
  m64030_190330_071939/102172020/ccs m64030_190330_071939/30016220/ccs -5450 99.8726 0 1217 10635 10635 1 946 10371 10371
  m64030_190330_071939/102172020/ccs m64030_190330_071939/28901470/ccs -3523 99.8517 0 4566 10635 10635 1 3541 9616 9723
  m64030_190330_071939/102172020/ccs m64030_190330_071939/49610760/ccs -2461 99.6697 0 0 4547 10635 1 0 4542 9121
  m64030_190330_071939/102172020/ccs m64030_190330_071939/52102340/ccs -1965 99.7537 0 6981 10635 10635 1 7165 10825 10825

Contained overlap.
  $ ${BIN_DIR}/raptor -x ovl-hifi -r ${PROJECT_DIR}/test-data/hifi-ovl/test-2-contained/reads.pile1.fasta -q ${PROJECT_DIR}/test-data/hifi-ovl/test-2-contained/reads.pile1.fasta -k 30 -w 80 -v 0 -t 1 -n 1 --out-fmt m4
  m64030_190330_071939/100008040/ccs m64030_190330_071939/117113640/ccs -5730 99.8781 0 0 9847 9847 0 388 10231 11163

Internal overlap.
  $ ${BIN_DIR}/raptor -x ovl-hifi -r ${PROJECT_DIR}/test-data/hifi-ovl/test-3-internal/reads.pile2.fasta -q ${PROJECT_DIR}/test-data/hifi-ovl/test-3-internal/reads.pile2.fasta -k 30 -w 80 -v 0 -t 1 -n 1 --out-fmt m4
  m64030_190330_071939/101451080/ccs m64030_190330_071939/161024730/ccs -2633 96.7319 0 2256 7060 8022 1 126 4930 12064
  m64030_190330_071939/101451080/ccs m64030_190330_071939/164825240/ccs -2505 96.2606 0 2227 7415 8022 0 5753 10965 11350
  m64030_190330_071939/101451080/ccs m64030_190330_071939/164823320/ccs -2314 95.9715 0 2256 7444 8022 1 1428 6637 10056
  m64030_190330_071939/101451080/ccs m64030_190330_071939/165414370/ccs -1392 95.4594 0 4513 7390 8022 0 0 2797 10706
  m64030_190330_071939/101451080/ccs m64030_190330_071939/160039920/ccs -500 98.9011 0 6261 7357 8022 1 7586 8678 8678