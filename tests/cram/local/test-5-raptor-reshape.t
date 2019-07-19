This runs raptor-reshape on a small input.
  $ ${BIN_DIR}/raptor-reshape -i ${PROJECT_DIR}/test-data/ecoli-small/reads.6x.fwd.fasta -o out --block-size 0.2 --split-blocks -v 0 && cat out.rdb
  V	0.2
  F	0	out.block.0	fasta
  S	0	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3005/0_5852	5852	0	0	5929
  S	1	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3414/0_11983	11983	0	5929	12061
  S	2	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292	24292	0	17990	24370
  S	3	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3981/0_5105	5105	0	42360	5182
  S	4	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4028/0_19001	19001	0	47542	19079
  S	5	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4477/0_7032	7032	0	66621	7109
  S	6	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4630/0_17021	17021	0	73730	17099
  S	7	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4630/17069_24358	7289	0	90829	7371
  S	8	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4722/0_15634	15634	0	98200	15712
  S	9	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4950/0_7575	7575	0	113912	7652
  S	10	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4980/0_20648	20648	0	121564	20726
  S	11	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/5314/0_12583	12583	0	142290	12661
  S	12	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/5838/0_5509	5509	0	154951	5586
  S	13	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/6299/0_5972	5972	0	160537	6049
  S	14	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/6323/0_15585	15585	0	166586	15663
  S	15	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/6651/0_6723	6723	0	182249	6800
  S	16	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/6708/0_6277	6277	0	189049	6354
  S	17	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/6867/6090_8949	2859	0	195403	2939
  S	18	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/6873/0_11335	11335	0	198342	11413
  S	19	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/7052/0_21860	21860	0	209755	21938
  B	0	0	20	230135
  F	1	out.block.1	fasta
  S	20	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/7060/15041_31664	16623	1	0	16705
  S	21	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/7239/0_23877	23877	1	16705	23955
  S	22	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/7239/23920_31197	7277	1	40660	7359
  S	23	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/7252/0_19635	19635	1	48019	19713
  S	24	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/7626/0_19169	19169	1	67732	19247
  S	25	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/7626/19215_23470	4255	1	86979	4337
  S	26	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/7642/0_6981	6981	1	91316	7058
  S	27	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/7768/0_20664	20664	1	98374	20742
  S	28	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/7768/20714_23917	3203	1	119116	3285
  S	29	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/7998/0_18132	18132	1	122401	18210
  S	30	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/8000/0_12334	12334	1	140611	12412
  S	31	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/8002/0_15731	15731	1	153023	15809
  S	32	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/8133/739_8981	8242	1	168832	8321
  S	33	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/8156/0_7011	7011	1	177153	7088
  S	34	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/8426/0_5361	5361	1	184241	5438
  S	35	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/8516/0_1610	1610	1	189679	1687
  S	36	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/8927/0_18304	18304	1	191366	18382
  S	37	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/8927/18346_35482	17136	1	209748	17218
  B	1	20	38	225545
  F	2	out.block.2	fasta
  S	38	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/8964/0_12026	12026	2	0	12104
  S	39	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/8976/3577_6456	2879	2	12104	2959
  S	40	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9104/0_17379	17379	2	15063	17457
  S	41	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9158/3746_18052	14306	2	32520	14387
  S	42	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9164/1351_15175	13824	2	46907	13905
  S	43	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9192/0_23164	23164	2	60812	23242
  S	44	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9301/0_23092	23092	2	84054	23170
  S	45	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9491/18951_24483	5532	2	107224	5614
  S	46	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9569/0_3326	3326	2	112838	3403
  S	47	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9569/3373_6916	3543	2	116241	3623
  S	48	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9569/6964_10108	3144	2	119864	3225
  S	49	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9569/10153_13253	3100	2	123089	3182
  S	50	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9569/13295_16423	3128	2	126271	3210
  S	51	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9569/16467_19630	3163	2	129481	3245
  S	52	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9569/19677_22845	3168	2	132726	3250
  S	53	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9569/22890_25973	3083	2	135976	3165
  S	54	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9569/26019_29177	3158	2	139141	3240
  S	55	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9569/29223_32322	3099	2	142381	3181
  S	56	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9569/32367_33943	1576	2	145562	1658
  S	57	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9767/9303_28629	19326	2	147220	19407
  S	58	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9916/0_15176	15176	2	166627	15254
  S	59	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9916/15222_18611	3389	2	181881	3471
  S	60	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/10301/0_5437	5437	2	185352	5515
  S	61	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/10307/0_4189	4189	2	190867	4267
  B	2	38	62	193207

Another test to check the edge case when one sequence overflows the block size, the block size is 0.09 and one seq is 100kb in length.
  $ ${BIN_DIR}/raptor-reshape -i ${PROJECT_DIR}/test-data/ecoli-small/double-ecoli-0-100000.fasta -o out --block-size 0.09 -v 0 && cat out.rdb
  V	0.2
  F	0	out.block.0	fasta
  S	0	1-gi|545778205|gb|U00096.3|	100000	0	0	100089
  B	0	0	1	100000
  S	1	2-gi|545778205|gb|U00096.3|	100000	0	100089	100089
  B	1	1	2	100000

Another test to check the edge case when one sequence overflows the block size, the block size is 0.11 and one seq is 100kb in length.
  $ ${BIN_DIR}/raptor-reshape -i ${PROJECT_DIR}/test-data/ecoli-small/double-ecoli-0-100000.fasta -o out --block-size 0.11 -v 0 && cat out.rdb
  V	0.2
  F	0	out.block.0	fasta
  S	0	1-gi|545778205|gb|U00096.3|	100000	0	0	100089
  S	1	2-gi|545778205|gb|U00096.3|	100000	0	100089	100089
  B	0	0	2	200000

Running Raptor using the RaptorDB as the reference file, loading all blocks.
  $ ${BIN_DIR}/raptor-reshape -i ${PROJECT_DIR}/test-data/ecoli-small/double-ecoli-0-100000.fasta -o out --block-size 0.09 -v 0 && ${BIN_DIR}/raptor -r out.rdb -q ${PROJECT_DIR}/test-data/ecoli-small/single_long_read.fa --out-fmt m4 -v 0
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292 1-gi|545778205|gb|U00096.3| -12323 51.5287 0 317 24259 24292 0 24806 47832 100000
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292 2-gi|545778205|gb|U00096.3| -12323 51.5287 0 317 24259 24292 0 24806 47832 100000

Running Raptor using the RaptorDB as the reference file, loading only block 0.
  $ ${BIN_DIR}/raptor-reshape -i ${PROJECT_DIR}/test-data/ecoli-small/double-ecoli-0-100000.fasta -o out --block-size 0.09 -v 0 && ${BIN_DIR}/raptor -r out.rdb -q ${PROJECT_DIR}/test-data/ecoli-small/single_long_read.fa --out-fmt m4 -v 0 --rdb-block-ref 0
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292 1-gi|545778205|gb|U00096.3| -12323 51.5287 0 317 24259 24292 0 24806 47832 100000

Running Raptor using the RaptorDB as the reference file, loading only block 1.
  $ ${BIN_DIR}/raptor-reshape -i ${PROJECT_DIR}/test-data/ecoli-small/double-ecoli-0-100000.fasta -o out --block-size 0.09 -v 0 && ${BIN_DIR}/raptor -r out.rdb -q ${PROJECT_DIR}/test-data/ecoli-small/single_long_read.fa --out-fmt m4 -v 0 --rdb-block-ref 1
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292 2-gi|545778205|gb|U00096.3| -12323 51.5287 0 317 24259 24292 0 24806 47832 100000

Running Raptor using the RaptorDB as the reference file, loading block out of range.
  $ ${BIN_DIR}/raptor-reshape -i ${PROJECT_DIR}/test-data/ecoli-small/double-ecoli-0-100000.fasta -o out --block-size 0.09 -v 0 && (${BIN_DIR}/raptor -r out.rdb -q ${PROJECT_DIR}/test-data/ecoli-small/single_long_read.fa --out-fmt m4 -v 0 --rdb-block-ref 5 2>&1 | sed 's/.*#5: //')
  Unexpected value found! Specified block_id is larger than the number of blocks in RaptorDB. block_id = 5, rasf->db_blocks().size() = 2.

Running Raptor using the RaptorDB as the reference file, ref block < 0, so all blocks should be processed.
  $ ${BIN_DIR}/raptor-reshape -i ${PROJECT_DIR}/test-data/ecoli-small/double-ecoli-0-100000.fasta -o out --block-size 0.09 -v 0 && (${BIN_DIR}/raptor -r out.rdb -q ${PROJECT_DIR}/test-data/ecoli-small/single_long_read.fa --out-fmt m4 -v 0 --rdb-block-ref -1)
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292 1-gi|545778205|gb|U00096.3| -12323 51.5287 0 317 24259 24292 0 24806 47832 100000
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292 2-gi|545778205|gb|U00096.3| -12323 51.5287 0 317 24259 24292 0 24806 47832 100000


Running Raptor using the RaptorDB as the reference file, loading all blocks.
  $ ${BIN_DIR}/raptor-reshape -i ${PROJECT_DIR}/test-data/ecoli-small/double-ecoli-0-100000.fasta -o out --block-size 0.09 -v 0 && ${BIN_DIR}/raptor -r out.rdb -q out.rdb --out-fmt m4 -v 0
  1-gi|545778205|gb|U00096.3| 1-gi|545778205|gb|U00096.3| -100000 100.014 0 0 99986 100000 0 0 99986 100000
  1-gi|545778205|gb|U00096.3| 2-gi|545778205|gb|U00096.3| -100000 100.014 0 0 99986 100000 0 0 99986 100000
  2-gi|545778205|gb|U00096.3| 1-gi|545778205|gb|U00096.3| -100000 100.014 0 0 99986 100000 0 0 99986 100000
  2-gi|545778205|gb|U00096.3| 2-gi|545778205|gb|U00096.3| -100000 100.014 0 0 99986 100000 0 0 99986 100000

Running Raptor using the RaptorDB as the reference file, ref block 0 and query block 0.
  $ ${BIN_DIR}/raptor-reshape -i ${PROJECT_DIR}/test-data/ecoli-small/double-ecoli-0-100000.fasta -o out --block-size 0.09 -v 0 && ${BIN_DIR}/raptor -r out.rdb -q out.rdb --out-fmt m4 -v 0 --rdb-block-ref 0 --rdb-block-query 0
  1-gi|545778205|gb|U00096.3| 1-gi|545778205|gb|U00096.3| -100000 100.014 0 0 99986 100000 0 0 99986 100000

Running Raptor using the RaptorDB as the reference file, ref block 0 and query block 1.
  $ ${BIN_DIR}/raptor-reshape -i ${PROJECT_DIR}/test-data/ecoli-small/double-ecoli-0-100000.fasta -o out --block-size 0.09 -v 0 && ${BIN_DIR}/raptor -r out.rdb -q out.rdb --out-fmt m4 -v 0 --rdb-block-ref 0 --rdb-block-query 1
  2-gi|545778205|gb|U00096.3| 1-gi|545778205|gb|U00096.3| -100000 100.014 0 0 99986 100000 0 0 99986 100000

Running Raptor using the RaptorDB as the reference file, entire ref, and query block 1.
  $ ${BIN_DIR}/raptor-reshape -i ${PROJECT_DIR}/test-data/ecoli-small/double-ecoli-0-100000.fasta -o out --block-size 0.09 -v 0 && ${BIN_DIR}/raptor -r out.rdb -q out.rdb --out-fmt m4 -v 0 --rdb-block-query 1
  2-gi|545778205|gb|U00096.3| 1-gi|545778205|gb|U00096.3| -100000 100.014 0 0 99986 100000 0 0 99986 100000
  2-gi|545778205|gb|U00096.3| 2-gi|545778205|gb|U00096.3| -100000 100.014 0 0 99986 100000 0 0 99986 100000

Running Raptor using the RaptorDB as the reference file, entire ref, and query block < 0 so all blocks should be processed.
  $ ${BIN_DIR}/raptor-reshape -i ${PROJECT_DIR}/test-data/ecoli-small/double-ecoli-0-100000.fasta -o out --block-size 0.09 -v 0 && (${BIN_DIR}/raptor -r out.rdb -q out.rdb --out-fmt m4 -v 0 --rdb-block-query -1)
  1-gi|545778205|gb|U00096.3| 1-gi|545778205|gb|U00096.3| -100000 100.014 0 0 99986 100000 0 0 99986 100000
  1-gi|545778205|gb|U00096.3| 2-gi|545778205|gb|U00096.3| -100000 100.014 0 0 99986 100000 0 0 99986 100000
  2-gi|545778205|gb|U00096.3| 1-gi|545778205|gb|U00096.3| -100000 100.014 0 0 99986 100000 0 0 99986 100000
  2-gi|545778205|gb|U00096.3| 2-gi|545778205|gb|U00096.3| -100000 100.014 0 0 99986 100000 0 0 99986 100000

Running Raptor using the RaptorDB as the reference file, entire ref, and query block out of range.
  $ ${BIN_DIR}/raptor-reshape -i ${PROJECT_DIR}/test-data/ecoli-small/double-ecoli-0-100000.fasta -o out --block-size 0.09 -v 0 && (${BIN_DIR}/raptor -r out.rdb -q out.rdb --out-fmt m4 -v 0 --rdb-block-query 5 2>&1 | sed 's/.*#5: //')
  Unexpected value found! The start_block_id >= num_blocks. start_block_id = 5, num_blocks = 2. Exiting.



#############################
This runs raptor-reshape on a small input, but symlinking the file instead of writing the data to new files.
  $ ${BIN_DIR}/raptor-reshape -i ${PROJECT_DIR}/test-data/ecoli-small/reads.6x.fwd.fasta -o out --block-size 0.2 --split-blocks -v 0 --symlink 2>&1 | grep "Option"
  Option '--split-blocks' cannot be used in combination with '--symlink'.

This runs raptor-reshape on a small input, but symlinking the file instead of writing the data to new files.
  $ (${BIN_DIR}/raptor-reshape -i ${PROJECT_DIR}/test-data/ecoli-small/reads.6x.fwd.fasta -o out --block-size 0.2 -v 0 --symlink && cat out.rdb | sed "s#${PROJECT_DIR}##g")
  V	0.2
  F	0	/test-data/ecoli-small/reads.6x.fwd.fasta	fasta
  S	0	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3005/0_5852	5852	0	0	5929
  S	1	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3414/0_11983	11983	0	5929	12061
  S	2	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292	24292	0	17990	24370
  S	3	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3981/0_5105	5105	0	42360	5182
  S	4	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4028/0_19001	19001	0	47542	19079
  S	5	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4477/0_7032	7032	0	66621	7109
  S	6	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4630/0_17021	17021	0	73730	17099
  S	7	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4630/17069_24358	7289	0	90829	7371
  S	8	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4722/0_15634	15634	0	98200	15712
  S	9	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4950/0_7575	7575	0	113912	7652
  S	10	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4980/0_20648	20648	0	121564	20726
  S	11	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/5314/0_12583	12583	0	142290	12661
  S	12	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/5838/0_5509	5509	0	154951	5586
  S	13	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/6299/0_5972	5972	0	160537	6049
  S	14	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/6323/0_15585	15585	0	166586	15663
  S	15	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/6651/0_6723	6723	0	182249	6800
  S	16	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/6708/0_6277	6277	0	189049	6354
  S	17	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/6867/6090_8949	2859	0	195403	2939
  S	18	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/6873/0_11335	11335	0	198342	11413
  S	19	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/7052/0_21860	21860	0	209755	21938
  B	0	0	20	230135
  S	20	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/7060/15041_31664	16623	0	231693	16705
  S	21	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/7239/0_23877	23877	0	248398	23955
  S	22	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/7239/23920_31197	7277	0	272353	7359
  S	23	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/7252/0_19635	19635	0	279712	19713
  S	24	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/7626/0_19169	19169	0	299425	19247
  S	25	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/7626/19215_23470	4255	0	318672	4337
  S	26	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/7642/0_6981	6981	0	323009	7058
  S	27	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/7768/0_20664	20664	0	330067	20742
  S	28	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/7768/20714_23917	3203	0	350809	3285
  S	29	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/7998/0_18132	18132	0	354094	18210
  S	30	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/8000/0_12334	12334	0	372304	12412
  S	31	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/8002/0_15731	15731	0	384716	15809
  S	32	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/8133/739_8981	8242	0	400525	8321
  S	33	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/8156/0_7011	7011	0	408846	7088
  S	34	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/8426/0_5361	5361	0	415934	5438
  S	35	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/8516/0_1610	1610	0	421372	1687
  S	36	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/8927/0_18304	18304	0	423059	18382
  S	37	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/8927/18346_35482	17136	0	441441	17218
  B	1	20	38	225545
  S	38	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/8964/0_12026	12026	0	458659	12104
  S	39	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/8976/3577_6456	2879	0	470763	2959
  S	40	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9104/0_17379	17379	0	473722	17457
  S	41	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9158/3746_18052	14306	0	491179	14387
  S	42	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9164/1351_15175	13824	0	505566	13905
  S	43	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9192/0_23164	23164	0	519471	23242
  S	44	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9301/0_23092	23092	0	542713	23170
  S	45	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9491/18951_24483	5532	0	565883	5614
  S	46	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9569/0_3326	3326	0	571497	3403
  S	47	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9569/3373_6916	3543	0	574900	3623
  S	48	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9569/6964_10108	3144	0	578523	3225
  S	49	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9569/10153_13253	3100	0	581748	3182
  S	50	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9569/13295_16423	3128	0	584930	3210
  S	51	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9569/16467_19630	3163	0	588140	3245
  S	52	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9569/19677_22845	3168	0	591385	3250
  S	53	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9569/22890_25973	3083	0	594635	3165
  S	54	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9569/26019_29177	3158	0	597800	3240
  S	55	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9569/29223_32322	3099	0	601040	3181
  S	56	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9569/32367_33943	1576	0	604221	1658
  S	57	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9767/9303_28629	19326	0	605879	19407
  S	58	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9916/0_15176	15176	0	625286	15254
  S	59	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9916/15222_18611	3389	0	640540	3471
  S	60	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/10301/0_5437	5437	0	644011	5515
  S	61	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/10307/0_4189	4189	0	649526	4266
  B	2	38	62	193207

Test symlink on a simple file.
  $ (${BIN_DIR}/raptor-reshape -i ${PROJECT_DIR}/test-data/ecoli-small/double-ecoli-0-100000.fasta -o out --block-size 0.09 --symlink -v 0 && cat out.rdb | sed "s#${PROJECT_DIR}##g")
  V	0.2
  F	0	/test-data/ecoli-small/double-ecoli-0-100000.fasta	fasta
  S	0	1-gi|545778205|gb|U00096.3|	100000	0	0	100089
  B	0	0	1	100000
  S	1	2-gi|545778205|gb|U00096.3|	100000	0	100089	100088
  B	1	1	2	100000

Test symlink on multiple input files.
  $ (${BIN_DIR}/raptor-reshape -i ${PROJECT_DIR}/test-data/ecoli-small/double-ecoli-0-100000.fasta -i ${PROJECT_DIR}/test-data/ecoli-small/double-ecoli-0-100000.fasta -o out --block-size 0.09 --symlink -v 0 && cat out.rdb | sed "s#${PROJECT_DIR}##g")
  V	0.2
  F	0	/test-data/ecoli-small/double-ecoli-0-100000.fasta	fasta
  F	1	/test-data/ecoli-small/double-ecoli-0-100000.fasta	fasta
  S	0	1-gi|545778205|gb|U00096.3|	100000	0	0	100089
  B	0	0	1	100000
  S	1	2-gi|545778205|gb|U00096.3|	100000	0	100089	100088
  B	1	1	2	100000
  S	2	1-gi|545778205|gb|U00096.3|	100000	1	0	100089
  B	2	2	3	100000
  S	3	2-gi|545778205|gb|U00096.3|	100000	1	100089	100088
  B	3	3	4	100000

Test renaming of sequences.
  $ ${BIN_DIR}/raptor-reshape -i ${PROJECT_DIR}/test-data/ecoli-small/double-ecoli-0-100000.fasta --rename -o out --block-size 0.09 -v 0 && cat out.rdb
  V	0.2
  F	0	out.block.0	fasta
  S	0	000000000	100000	0	0	100012
  B	0	0	1	100000
  S	1	000000001	100000	0	100012	100012
  B	1	1	2	100000

Test the FOFN functionality.
  $ (echo "${PROJECT_DIR}/test-data/ecoli-small/double-ecoli-0-100000.fasta" > input.fofn) && (echo "${PROJECT_DIR}/test-data/ecoli-small/double-ecoli-0-100000.fasta" >> input.fofn) && (${BIN_DIR}/raptor-reshape -i input.fofn -o out --block-size 0.09 --symlink -v 0 && cat out.rdb | sed "s#${PROJECT_DIR}##g")
  V	0.2
  F	0	/test-data/ecoli-small/double-ecoli-0-100000.fasta	fasta
  F	1	/test-data/ecoli-small/double-ecoli-0-100000.fasta	fasta
  S	0	1-gi|545778205|gb|U00096.3|	100000	0	0	100089
  B	0	0	1	100000
  S	1	2-gi|545778205|gb|U00096.3|	100000	0	100089	100088
  B	1	1	2	100000
  S	2	1-gi|545778205|gb|U00096.3|	100000	1	0	100089
  B	2	2	3	100000
  S	3	2-gi|545778205|gb|U00096.3|	100000	1	100089	100088
  B	3	3	4	100000

Test the PacBio XML functionality.
  $ (ln -s ${PROJECT_DIR}/test-data/small-xml/subreads1.bam && ln -s ${PROJECT_DIR}/test-data/small-xml/subreads1.bam.pbi && ln -s ${PROJECT_DIR}/test-data/small-xml/subreads2.bam && ln -s ${PROJECT_DIR}/test-data/small-xml/subreads2.bam.pbi && ln -s ${PROJECT_DIR}/test-data/small-xml/subreads3.bam && ln -s ${PROJECT_DIR}/test-data/small-xml/subreads3.bam.pbi && ln -s ${PROJECT_DIR}/test-data/small-xml/subreadset.xml) && ${BIN_DIR}/raptor-reshape -i subreadset.xml --symlink -o out -v 0 && cat out.rdb
  V	0.2
  F	0	./subreads1.bam	bam
  F	1	./subreads2.bam	bam
  F	2	./subreads3.bam	bam
  S	0	ref1/1/0_42148	42148	0	29032448	865206272
  S	1	ref1/2/0_42124	42124	0	894238720	864354304
  S	2	ref1/3/0_42176	42176	0	1758593024	866451456
  S	3	ref1/4/0_42179	42179	0	2625044480	866123776
  S	4	ref1/5/0_42280	42280	0	3491168256	867565568
  S	5	ref2/1/0_42148	42148	1	29032448	865861632
  S	6	ref2/2/0_42124	42124	1	894894080	861798400
  S	7	ref2/3/0_42176	42176	1	1756692480	865206272
  S	8	ref2/4/0_42179	42179	1	2621898752	866320384
  S	9	ref2/5/0_42280	42280	1	3488219136	866451456
  S	10	ref1/158/0_23696	23696	2	29032448	490733568
  S	11	ref1/159/0_42173	42173	2	519766016	866451456
  S	12	ref1/160/0_42090	42090	2	1386217472	862978048
  S	13	ref1/161/0_38694	38694	2	2249195520	796917760
  S	14	ref1/162/0_42180	42180	2	3046113280	865468416
  B	0	0	15	610647
