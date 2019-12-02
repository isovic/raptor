Test overlapping of blocks with single arcs. The "-x ovl-raw" implicitly defines "--overlap-skip-self". Blocks 0 and 0 are compared.
  $ (${BIN_DIR}/raptor-reshape -v 0 -i ${PROJECT_DIR}/test-data/ecoli-small/two_identical_long_reads.fa --block-size 0.01 -o reads --symlink) &&  ${BIN_DIR}/raptor --out-fmt paf -r reads.rdb -q reads.rdb -v 0 -x ovl-raw -k 15 -w 5 --overlap-single-arc --rdb-block-ref 0 --rdb-block-query 0

Test overlapping of blocks with single arcs. The "-x ovl-raw" implicitly defines "--overlap-skip-self". Blocks 0 and 1 are compared.
  $ (${BIN_DIR}/raptor-reshape -v 0 -i ${PROJECT_DIR}/test-data/ecoli-small/two_identical_long_reads.fa --block-size 0.01 -o reads --symlink) &&  ${BIN_DIR}/raptor --out-fmt paf -r reads.rdb -q reads.rdb -v 0 -x ovl-raw -k 15 -w 5 --overlap-single-arc --rdb-block-ref 0 --rdb-block-query 1
  read-2	24292	2	24275	+	read-1	24292	2	24275	24268	24273	255	cm:i:1756	NM:i:-1	AS:i:24268	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*

Test overlapping of blocks with single arcs. The "-x ovl-raw" implicitly defines "--overlap-skip-self". Blocks 1 and 0 are compared.
  $ (${BIN_DIR}/raptor-reshape -v 0 -i ${PROJECT_DIR}/test-data/ecoli-small/two_identical_long_reads.fa --block-size 0.01 -o reads --symlink) &&  ${BIN_DIR}/raptor --out-fmt paf -r reads.rdb -q reads.rdb -v 0 -x ovl-raw -k 15 -w 5 --overlap-single-arc --rdb-block-ref 1 --rdb-block-query 0

Test overlapping of blocks with single arcs. The "-x ovl-raw" implicitly defines "--overlap-skip-self". Blocks 1 and 1 are compared.
  $ (${BIN_DIR}/raptor-reshape -v 0 -i ${PROJECT_DIR}/test-data/ecoli-small/two_identical_long_reads.fa --block-size 0.01 -o reads --symlink) &&  ${BIN_DIR}/raptor --out-fmt paf -r reads.rdb -q reads.rdb -v 0 -x ovl-raw -k 15 -w 5 --overlap-single-arc --rdb-block-ref 1 --rdb-block-query 1



Test overlapping of blocks with dual arcs. The "-x ovl-raw" implicitly defines "--overlap-skip-self". Blocks 0 and 0 are compared.
  $ (${BIN_DIR}/raptor-reshape -v 0 -i ${PROJECT_DIR}/test-data/ecoli-small/two_identical_long_reads.fa --block-size 0.01 -o reads --symlink) &&  ${BIN_DIR}/raptor --out-fmt paf -r reads.rdb -q reads.rdb -v 0 -x ovl-raw -k 15 -w 5 --rdb-block-ref 0 --rdb-block-query 0

Test overlapping of blocks with dual arcs. The "-x ovl-raw" implicitly defines "--overlap-skip-self". Blocks 0 and 1 are compared.
  $ (${BIN_DIR}/raptor-reshape -v 0 -i ${PROJECT_DIR}/test-data/ecoli-small/two_identical_long_reads.fa --block-size 0.01 -o reads --symlink) &&  ${BIN_DIR}/raptor --out-fmt paf -r reads.rdb -q reads.rdb -v 0 -x ovl-raw -k 15 -w 5 --rdb-block-ref 0 --rdb-block-query 1
  read-2	24292	2	24275	+	read-1	24292	2	24275	24268	24273	255	cm:i:1756	NM:i:-1	AS:i:24268	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*

Test overlapping of blocks with dual arcs. The "-x ovl-raw" implicitly defines "--overlap-skip-self". Blocks 1 and 0 are compared.
  $ (${BIN_DIR}/raptor-reshape -v 0 -i ${PROJECT_DIR}/test-data/ecoli-small/two_identical_long_reads.fa --block-size 0.01 -o reads --symlink) &&  ${BIN_DIR}/raptor --out-fmt paf -r reads.rdb -q reads.rdb -v 0 -x ovl-raw -k 15 -w 5 --rdb-block-ref 1 --rdb-block-query 0
  read-1	24292	2	24275	+	read-2	24292	2	24275	24268	24273	255	cm:i:1756	NM:i:-1	AS:i:24268	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*

Test overlapping of blocks with dual arcs. The "-x ovl-raw" implicitly defines "--overlap-skip-self". Blocks 1 and 1 are compared.
  $ (${BIN_DIR}/raptor-reshape -v 0 -i ${PROJECT_DIR}/test-data/ecoli-small/two_identical_long_reads.fa --block-size 0.01 -o reads --symlink) &&  ${BIN_DIR}/raptor --out-fmt paf -r reads.rdb -q reads.rdb -v 0 -x ovl-raw -k 15 -w 5 --rdb-block-ref 1 --rdb-block-query 1



Test overlapping on PacBio XML input.
  $ ${BIN_DIR}/raptor -x ovl-raw -r ${PROJECT_DIR}/test-data/small-xml/subreadset.xml -q ${PROJECT_DIR}/test-data/small-xml/subreadset.xml -v 0 | cut -f1-11
  ref1/1/0_42148	42148	624	41200	+	ref2/1/0_42148	42148	1408	41935	4685	40527
  ref1/1/0_42148	42148	20186	41924	+	ref1/5/0_42280	42280	121	22040	2912	21919
  ref1/1/0_42148	42148	15041	42067	+	ref2/5/0_42280	42280	168	27401	2754	27233
  ref1/1/0_42148	42148	35164	42080	+	ref1/162/0_42180	42180	173	7228	952	7055
  ref1/2/0_42124	42124	3439	42049	+	ref2/2/0_42124	42124	457	39152	5136	38695
  ref1/2/0_42124	42124	157	30009	+	ref2/4/0_42179	42179	12317	42110	3745	29793
  ref1/2/0_42124	42124	69	25890	+	ref1/4/0_42179	42179	16372	42096	3589	25724
  ref1/2/0_42124	42124	157	28945	+	ref2/3/0_42176	42176	13267	42043	3573	28776
  ref1/2/0_42124	42124	99	24780	+	ref1/3/0_42176	42176	17451	42063	3367	24612
  ref1/2/0_42124	42124	674	18130	+	ref1/159/0_42173	42173	24835	42151	2528	17316
  ref1/3/0_42176	42176	1489	42149	+	ref1/4/0_42179	42179	423	41086	5321	40663
  ref1/3/0_42176	42176	208	35461	+	ref1/159/0_42173	42173	7004	42151	4655	35147
  ref1/3/0_42176	42176	11177	42118	+	ref2/3/0_42176	42176	6961	37935	3754	30974
  ref1/3/0_42176	42176	17451	42063	+	ref1/2/0_42124	42124	99	24780	3367	24681
  ref1/3/0_42176	42176	11176	42118	+	ref2/4/0_42179	42179	6019	36957	3152	30938
  ref1/3/0_42176	42176	20490	42101	+	ref2/2/0_42124	42124	157	21841	2647	21684
  ref1/4/0_42179	42179	423	41086	+	ref1/3/0_42176	42176	1489	42149	5321	40660
  ref1/4/0_42179	42179	41	34345	+	ref1/159/0_42173	42173	7910	42096	4364	34186
  ref1/4/0_42179	42179	10462	42089	+	ref2/4/0_42179	42179	6290	38003	3701	31713
  ref1/4/0_42179	42179	16372	42096	+	ref1/2/0_42124	42124	69	25890	3590	25821
  ref1/4/0_42179	42179	10351	42096	+	ref2/3/0_42176	42176	7120	38985	3529	31865
  ref1/4/0_42179	42179	19591	42008	+	ref2/2/0_42124	42124	317	22827	2745	22510
  ref1/5/0_42280	42280	94	37049	+	ref2/5/0_42280	42280	5362	42208	3801	36846
  ref1/5/0_42280	42280	15148	42255	+	ref1/162/0_42180	42180	157	27359	3558	27202
  ref1/5/0_42280	42280	121	22040	+	ref1/1/0_42148	42148	20186	41924	2912	21738
  ref1/5/0_42280	42280	442	21410	+	ref2/1/0_42148	42148	21270	42040	2098	20770
  ref1/5/0_42280	42280	33300	42139	+	ref1/160/0_42090	42090	34	8907	1188	8873
  ref2/1/0_42148	42148	1408	41935	+	ref1/1/0_42148	42148	624	41200	4684	40576
  ref2/1/0_42148	42148	15785	42062	+	ref2/5/0_42280	42280	168	26656	3786	26488
  ref2/1/0_42148	42148	21270	42040	+	ref1/5/0_42280	42280	442	21410	2098	20968
  ref2/1/0_42148	42148	36077	42040	+	ref1/162/0_42180	42180	329	6434	819	6105
  ref2/2/0_42124	42124	457	39152	+	ref1/2/0_42124	42124	3439	42049	5136	38610
  ref2/2/0_42124	42124	157	27045	+	ref2/4/0_42179	42179	15290	42097	3254	26807
  ref2/2/0_42124	42124	512	26108	+	ref2/3/0_42176	42176	16619	42145	3033	25526
  ref2/2/0_42124	42124	317	22827	+	ref1/4/0_42179	42179	19591	42008	2747	22417
  ref2/2/0_42124	42124	157	21841	+	ref1/3/0_42176	42176	20490	42101	2647	21611
  ref2/2/0_42124	42124	121	15107	+	ref1/159/0_42173	42173	27230	42114	2142	14884
  ref2/3/0_42176	42176	1505	42145	+	ref2/4/0_42179	42179	516	41155	4911	40639
  ref2/3/0_42176	42176	6961	37935	+	ref1/3/0_42176	42176	11177	42118	3754	30941
  ref2/3/0_42176	42176	13267	42043	+	ref1/2/0_42124	42124	157	28945	3573	28788
  ref2/3/0_42176	42176	7120	38985	+	ref1/4/0_42179	42179	10351	42096	3528	31745
  ref2/3/0_42176	42176	16619	42145	+	ref2/2/0_42124	42124	512	26108	3031	25596
  ref2/3/0_42176	42176	6962	31277	+	ref1/159/0_42173	42173	18019	42142	2884	24123
  ref2/3/0_42176	42176	209	2209	+	ref1/160/0_42090	42090	40049	42047	295	1998
  ref2/4/0_42179	42179	516	41155	+	ref2/3/0_42176	42176	1505	42145	4912	40640
  ref2/4/0_42179	42179	12317	42110	+	ref1/2/0_42124	42124	157	30009	3745	29852
  ref2/4/0_42179	42179	6290	38003	+	ref1/4/0_42179	42179	10462	42089	3701	31627
  ref2/4/0_42179	42179	6021	30150	+	ref1/159/0_42173	42173	18019	42059	3381	24040
  ref2/4/0_42179	42179	15290	42097	+	ref2/2/0_42124	42124	157	27045	3252	26888
  ref2/4/0_42179	42179	6019	36957	+	ref1/3/0_42176	42176	11176	42118	3152	30942
  ref2/4/0_42179	42179	36	1213	+	ref1/160/0_42090	42090	40874	42021	161	1147
  ref2/5/0_42280	42280	5362	42208	+	ref1/5/0_42280	42280	94	37049	3801	36955
  ref2/5/0_42280	42280	168	26656	+	ref2/1/0_42148	42148	15785	42062	3785	26277
  ref2/5/0_42280	42280	20395	42106	+	ref1/162/0_42180	42180	157	21990	2752	21833
  ref2/5/0_42280	42280	168	27401	+	ref1/1/0_42148	42148	15041	42067	2752	27026
  ref2/5/0_42280	42280	38780	41994	+	ref1/160/0_42090	42090	333	3556	387	3223
  ref1/158/0_23696	23696	169	23441	+	ref1/161/0_38694	38694	15161	38443	2705	23282
  ref1/159/0_42173	42173	7004	42151	+	ref1/3/0_42176	42176	208	35461	4654	35253
  ref1/159/0_42173	42173	7910	42096	+	ref1/4/0_42179	42179	41	34345	4363	34304
  ref1/159/0_42173	42173	18019	42059	+	ref2/4/0_42179	42179	6021	30150	3381	24129
  ref1/159/0_42173	42173	18019	42142	+	ref2/3/0_42176	42176	6962	31277	2884	24315
  ref1/159/0_42173	42173	24835	42151	+	ref1/2/0_42124	42124	674	18130	2528	17456
  ref1/159/0_42173	42173	27230	42114	+	ref2/2/0_42124	42124	121	15107	2142	14986
  ref1/160/0_42090	42090	40	23910	+	ref1/162/0_42180	42180	18372	42140	3080	23768
  ref1/160/0_42090	42090	34	8907	+	ref1/5/0_42280	42280	33300	42139	1188	8839
  ref1/160/0_42090	42090	333	3556	+	ref2/5/0_42280	42280	38780	41994	387	3214
  ref1/160/0_42090	42090	40049	42047	+	ref2/3/0_42176	42176	209	2209	295	2000
  ref1/160/0_42090	42090	40874	42021	+	ref2/4/0_42179	42179	36	1213	161	1177
  ref1/161/0_38694	38694	15161	38443	+	ref1/158/0_23696	23696	169	23441	2705	23272
  ref1/162/0_42180	42180	157	27359	+	ref1/5/0_42280	42280	15148	42255	3559	27107
  ref1/162/0_42180	42180	18372	42140	+	ref1/160/0_42090	42090	40	23910	3080	23870
  ref1/162/0_42180	42180	157	21990	+	ref2/5/0_42280	42280	20395	42106	2752	21711
  ref1/162/0_42180	42180	173	7228	+	ref1/1/0_42148	42148	35164	42080	952	6916
  ref1/162/0_42180	42180	329	6434	+	ref2/1/0_42148	42148	36077	42040	819	5963
