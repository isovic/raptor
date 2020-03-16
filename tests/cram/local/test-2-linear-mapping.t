Mappping the same pair of reads where one is query and the other target, and then reversing the pair (i.e. overlapping with symmetric arcs). This should return identical results.
  $ ${BIN_DIR}/raptor -x ovl-raw -k 15 -w 5 --min-map-len 100 --out-fmt m4 -r ${PROJECT_DIR}/test-data/consistent-query-target-results/reads.14183.and.16026.fasta -q ${PROJECT_DIR}/test-data/consistent-query-target-results/reads.14183.and.16026.fasta -v 0
  14183 16026 -842 100.0000 0 44 884 19004 0 8831 9674 9689
  16026 14183 -842 100.0000 0 8831 9674 9689 0 44 884 19004

Mapping first 100 reads, all of which should be in fwd.
  $ ${BIN_DIR}/raptor -n 100 -r ${PROJECT_DIR}/test-data/ecoli-small/ecoli-0-100000.fasta -q ${PROJECT_DIR}/test-data/ecoli-small/reads.6x.fwd.fasta -v 0 | cut -f1-11,13-
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3005/0_5852	5852	19	5705	+	gi|545778205|gb|U00096.3|	100000	40536	46052	3201	5516	cm:i:274	NM:i:-1	AS:i:3195	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3414/0_11983	11983	10949	11786	+	gi|545778205|gb|U00096.3|	100000	19806	20537	239	731	cm:i:21	NM:i:-1	AS:i:239	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292	24292	317	24259	+	gi|545778205|gb|U00096.3|	100000	24806	47832	12337	23026	cm:i:1050	NM:i:-1	AS:i:12323	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3981/0_5105	5105	30	4716	+	gi|545778205|gb|U00096.3|	100000	88252	92731	1969	4479	cm:i:175	NM:i:-1	AS:i:1965	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4028/0_19001	19001	5897	6375	+	gi|545778205|gb|U00096.3|	100000	19973	20417	73	444	cm:i:6	NM:i:-1	AS:i:73	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4477/0_7032	7032	1302	2054	+	gi|545778205|gb|U00096.3|	100000	19806	20537	396	731	cm:i:33	NM:i:-1	AS:i:396	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4630/0_17021	17021	6257	6825	+	gi|545778205|gb|U00096.3|	100000	19998	20544	185	546	cm:i:17	NM:i:-1	AS:i:185	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4630/17069_24358	7289	6060	6802	+	gi|545778205|gb|U00096.3|	100000	19816	20532	85	716	cm:i:8	NM:i:-1	AS:i:85	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4722/0_15634	15634	36	15605	+	gi|545778205|gb|U00096.3|	100000	41228	56093	8771	14865	cm:i:751	NM:i:-1	AS:i:8765	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4950/0_7575	7575	397	705	+	gi|545778205|gb|U00096.3|	100000	16382	16691	80	309	cm:i:7	NM:i:-1	AS:i:80	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4980/0_20648	20648	1788	2561	+	gi|545778205|gb|U00096.3|	100000	19810	20537	396	727	cm:i:36	NM:i:-1	AS:i:396	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/5314/0_12583	12583	4606	5322	+	gi|545778205|gb|U00096.3|	100000	19831	20550	371	719	cm:i:35	NM:i:-1	AS:i:371	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/5838/0_5509	5509	1873	5491	+	gi|545778205|gb|U00096.3|	100000	55535	58780	1338	3245	cm:i:120	NM:i:-1	AS:i:1337	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/6299/0_5972	5972	201	5877	+	gi|545778205|gb|U00096.3|	100000	44452	49837	2736	5385	cm:i:238	NM:i:-1	AS:i:2736	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/6323/0_15585	15585	11858	12538	+	gi|545778205|gb|U00096.3|	100000	19857	20529	385	672	cm:i:32	NM:i:-1	AS:i:385	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/6651/0_6723	6723	6	6374	+	gi|545778205|gb|U00096.3|	100000	90809	97181	1615	6372	cm:i:142	NM:i:-1	AS:i:1615	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/6708/0_6277	6277	50	5932	+	gi|545778205|gb|U00096.3|	100000	18928	24576	3180	5648	cm:i:276	NM:i:-1	AS:i:3176	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/6867/6090_8949	2859	2525	2819	+	gi|545778205|gb|U00096.3|	100000	19802	20096	184	294	cm:i:15	NM:i:-1	AS:i:184	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/6873/0_11335	11335	6589	7338	+	gi|545778205|gb|U00096.3|	100000	19831	20537	317	706	cm:i:27	NM:i:-1	AS:i:315	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/7052/0_21860	21860	15720	21410	+	gi|545778205|gb|U00096.3|	100000	64	5577	2989	5513	cm:i:262	NM:i:-1	AS:i:2988	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/7060/15041_31664	16623	6948	8439	+	gi|545778205|gb|U00096.3|	100000	15435	16714	573	1279	cm:i:49	NM:i:-1	AS:i:573	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/7239/0_23877	23877	5016	6190	+	gi|545778205|gb|U00096.3|	100000	15441	16623	608	1182	cm:i:55	NM:i:-1	AS:i:608	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/7239/23920_31197	7277	5088	6275	+	gi|545778205|gb|U00096.3|	100000	15493	16699	728	1206	cm:i:58	NM:i:-1	AS:i:728	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/7252/0_19635	19635	1977	3263	+	gi|545778205|gb|U00096.3|	100000	15389	16714	552	1325	cm:i:50	NM:i:-1	AS:i:552	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/7626/0_19169	19169	397	19150	+	gi|545778205|gb|U00096.3|	100000	80037	98362	8576	18325	cm:i:745	NM:i:-1	AS:i:8576	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/7626/19215_23470	4255	30	4226	+	gi|545778205|gb|U00096.3|	100000	94222	98351	1979	4129	cm:i:171	NM:i:-1	AS:i:1977	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/7642/0_6981	6981	3444	6664	+	gi|545778205|gb|U00096.3|	100000	119	3443	673	3324	cm:i:61	NM:i:-1	AS:i:673	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/7768/0_20664	20664	187	20650	+	gi|545778205|gb|U00096.3|	100000	78542	97924	10359	19382	cm:i:877	NM:i:-1	AS:i:10354	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/7768/20714_23917	3203	48	3187	+	gi|545778205|gb|U00096.3|	100000	94911	97924	1508	3013	cm:i:130	NM:i:-1	AS:i:1507	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/7998/0_18132	18132	486	1138	+	gi|545778205|gb|U00096.3|	100000	19797	20534	88	737	cm:i:9	NM:i:-1	AS:i:88	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/8000/0_12334	12334	275	12289	+	gi|545778205|gb|U00096.3|	100000	30187	42417	2981	12230	cm:i:267	NM:i:-1	AS:i:2981	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/8002/0_15731	15731	15017	15633	+	gi|545778205|gb|U00096.3|	100000	19796	20332	202	536	cm:i:19	NM:i:-1	AS:i:202	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/8133/739_8981	8242	26	8228	+	gi|545778205|gb|U00096.3|	100000	14584	22488	3911	7904	cm:i:342	NM:i:-1	AS:i:3904	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/8156/0_7011	7011	62	6913	+	gi|545778205|gb|U00096.3|	100000	1290	7544	3075	6254	cm:i:259	NM:i:-1	AS:i:3073	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/8426/0_5361	5361	4251	5341	+	gi|545778205|gb|U00096.3|	100000	15399	16493	514	1094	cm:i:47	NM:i:-1	AS:i:514	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/8516/0_1610	1610	47	1230	+	gi|545778205|gb|U00096.3|	100000	82445	83579	341	1134	cm:i:27	NM:i:-1	AS:i:341	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/8927/0_18304	18304	4992	6288	+	gi|545778205|gb|U00096.3|	100000	15404	16699	638	1295	cm:i:57	NM:i:-1	AS:i:638	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/8927/18346_35482	17136	4070	5371	+	gi|545778205|gb|U00096.3|	100000	15419	16694	457	1275	cm:i:42	NM:i:-1	AS:i:457	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/8964/0_12026	12026	11413	12011	+	gi|545778205|gb|U00096.3|	100000	19826	20406	173	580	cm:i:15	NM:i:-1	AS:i:171	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/8976/3577_6456	2879	14	2816	+	gi|545778205|gb|U00096.3|	100000	52591	55442	1013	2851	cm:i:87	NM:i:-1	AS:i:1013	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9104/0_17379	17379	1316	1997	+	gi|545778205|gb|U00096.3|	100000	19836	20537	277	701	cm:i:26	NM:i:-1	AS:i:277	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9158/3746_18052	14306	12	14111	+	gi|545778205|gb|U00096.3|	100000	57994	71219	4001	13225	cm:i:356	NM:i:-1	AS:i:3998	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9164/1351_15175	13824	8853	9628	+	gi|545778205|gb|U00096.3|	100000	19796	20547	450	751	cm:i:37	NM:i:-1	AS:i:449	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9192/0_23164	23164	11344	12100	+	gi|545778205|gb|U00096.3|	100000	19826	20537	381	711	cm:i:31	NM:i:-1	AS:i:381	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9301/0_23092	23092	20283	21047	+	gi|545778205|gb|U00096.3|	100000	19797	20532	468	735	cm:i:41	NM:i:-1	AS:i:468	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9491/18951_24483	5532	8	475	+	gi|545778205|gb|U00096.3|	100000	99411	99851	195	440	cm:i:17	NM:i:-1	AS:i:195	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9569/0_3326	3326	38	3033	+	gi|545778205|gb|U00096.3|	100000	40503	43339	1303	2836	cm:i:115	NM:i:-1	AS:i:1300	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9569/3373_6916	3543	617	3527	+	gi|545778205|gb|U00096.3|	100000	40853	43498	1338	2645	cm:i:114	NM:i:-1	AS:i:1337	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9569/6964_10108	3144	43	3112	+	gi|545778205|gb|U00096.3|	100000	40534	43480	1717	2946	cm:i:148	NM:i:-1	AS:i:1717	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9569/10153_13253	3100	44	3047	+	gi|545778205|gb|U00096.3|	100000	40536	43465	1700	2929	cm:i:142	NM:i:-1	AS:i:1700	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9569/13295_16423	3128	34	3088	+	gi|545778205|gb|U00096.3|	100000	40523	43478	1649	2955	cm:i:138	NM:i:-1	AS:i:1649	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9569/16467_19630	3163	4	3146	+	gi|545778205|gb|U00096.3|	100000	40497	43497	1683	3000	cm:i:140	NM:i:-1	AS:i:1683	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9569/19677_22845	3168	4	3151	+	gi|545778205|gb|U00096.3|	100000	40497	43497	1714	3000	cm:i:147	NM:i:-1	AS:i:1710	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9569/22890_25973	3083	24	3039	+	gi|545778205|gb|U00096.3|	100000	40515	43474	1700	2959	cm:i:149	NM:i:-1	AS:i:1700	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9569/26019_29177	3158	68	3002	+	gi|545778205|gb|U00096.3|	100000	40555	43369	1366	2814	cm:i:118	NM:i:-1	AS:i:1365	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9569/29223_32322	3099	46	3082	+	gi|545778205|gb|U00096.3|	100000	40534	43497	1730	2963	cm:i:150	NM:i:-1	AS:i:1728	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9569/32367_33943	1576	14	1549	+	gi|545778205|gb|U00096.3|	100000	42049	43489	774	1440	cm:i:68	NM:i:-1	AS:i:773	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9767/9303_28629	19326	24	19236	+	gi|545778205|gb|U00096.3|	100000	41488	60123	7588	18635	cm:i:657	NM:i:-1	AS:i:7588	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9916/0_15176	15176	1145	14835	+	gi|545778205|gb|U00096.3|	100000	3	12918	6538	12915	cm:i:570	NM:i:-1	AS:i:6533	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9916/15222_18611	3389	1090	3364	+	gi|545778205|gb|U00096.3|	100000	3	2243	1147	2240	cm:i:101	NM:i:-1	AS:i:1147	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/10301/0_5437	5437	0	5255	+	gi|545778205|gb|U00096.3|	100000	90091	95196	2491	5105	cm:i:215	NM:i:-1	AS:i:2490	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/10307/0_4189	4189	53	3813	+	gi|545778205|gb|U00096.3|	100000	3605	7155	1841	3550	cm:i:155	NM:i:-1	AS:i:1838	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*

Mapping first 100 reads, all of which should be in rev.
  $ ${BIN_DIR}/raptor -n 100 -r ${PROJECT_DIR}/test-data/ecoli-small/ecoli-0-100000.fasta -q ${PROJECT_DIR}/test-data/ecoli-small/reads.6x.rev.fasta -v 0 | cut -f1-11,13-
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3005/0_5852	5852	147	5833	-	gi|545778205|gb|U00096.3|	100000	40536	46052	3199	5516	cm:i:274	NM:i:-1	AS:i:3195	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3414/0_11983	11983	197	1034	-	gi|545778205|gb|U00096.3|	100000	19806	20537	239	731	cm:i:21	NM:i:-1	AS:i:239	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292	24292	33	23975	-	gi|545778205|gb|U00096.3|	100000	24806	47832	12337	23026	cm:i:1050	NM:i:-1	AS:i:12324	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3981/0_5105	5105	389	5075	-	gi|545778205|gb|U00096.3|	100000	88252	92731	1968	4479	cm:i:175	NM:i:-1	AS:i:1965	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4028/0_19001	19001	12626	13104	-	gi|545778205|gb|U00096.3|	100000	19973	20417	73	444	cm:i:6	NM:i:-1	AS:i:73	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4477/0_7032	7032	4978	5730	-	gi|545778205|gb|U00096.3|	100000	19806	20537	396	731	cm:i:33	NM:i:-1	AS:i:396	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4630/0_17021	17021	10196	10764	-	gi|545778205|gb|U00096.3|	100000	19998	20544	185	546	cm:i:17	NM:i:-1	AS:i:185	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4630/17069_24358	7289	487	1229	-	gi|545778205|gb|U00096.3|	100000	19816	20532	85	716	cm:i:8	NM:i:-1	AS:i:85	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4722/0_15634	15634	29	15598	-	gi|545778205|gb|U00096.3|	100000	41228	56093	8776	14865	cm:i:751	NM:i:-1	AS:i:8766	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4950/0_7575	7575	6870	7178	-	gi|545778205|gb|U00096.3|	100000	16382	16691	80	309	cm:i:7	NM:i:-1	AS:i:80	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4980/0_20648	20648	18087	18860	-	gi|545778205|gb|U00096.3|	100000	19810	20537	396	727	cm:i:36	NM:i:-1	AS:i:396	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/5314/0_12583	12583	7261	7977	-	gi|545778205|gb|U00096.3|	100000	19831	20550	371	719	cm:i:35	NM:i:-1	AS:i:371	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/5838/0_5509	5509	18	3636	-	gi|545778205|gb|U00096.3|	100000	55535	58780	1339	3245	cm:i:120	NM:i:-1	AS:i:1338	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/6299/0_5972	5972	95	5771	-	gi|545778205|gb|U00096.3|	100000	44452	49837	2737	5385	cm:i:238	NM:i:-1	AS:i:2736	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/6323/0_15585	15585	3047	3727	-	gi|545778205|gb|U00096.3|	100000	19857	20529	385	672	cm:i:32	NM:i:-1	AS:i:385	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/6651/0_6723	6723	349	6717	-	gi|545778205|gb|U00096.3|	100000	90809	97181	1615	6372	cm:i:142	NM:i:-1	AS:i:1615	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/6708/0_6277	6277	345	6227	-	gi|545778205|gb|U00096.3|	100000	18928	24576	3181	5648	cm:i:276	NM:i:-1	AS:i:3176	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/6867/6090_8949	2859	40	334	-	gi|545778205|gb|U00096.3|	100000	19802	20096	184	294	cm:i:15	NM:i:-1	AS:i:184	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/6873/0_11335	11335	3997	4746	-	gi|545778205|gb|U00096.3|	100000	19831	20537	317	706	cm:i:27	NM:i:-1	AS:i:315	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/7052/0_21860	21860	450	6140	-	gi|545778205|gb|U00096.3|	100000	64	5577	2990	5513	cm:i:262	NM:i:-1	AS:i:2990	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/7060/15041_31664	16623	8184	9675	-	gi|545778205|gb|U00096.3|	100000	15435	16714	573	1279	cm:i:49	NM:i:-1	AS:i:573	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/7239/0_23877	23877	17687	18861	-	gi|545778205|gb|U00096.3|	100000	15441	16623	608	1182	cm:i:55	NM:i:-1	AS:i:608	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/7239/23920_31197	7277	1002	2189	-	gi|545778205|gb|U00096.3|	100000	15493	16699	728	1206	cm:i:58	NM:i:-1	AS:i:728	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/7252/0_19635	19635	16372	17658	-	gi|545778205|gb|U00096.3|	100000	15389	16714	552	1325	cm:i:50	NM:i:-1	AS:i:552	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/7626/0_19169	19169	19	18772	-	gi|545778205|gb|U00096.3|	100000	80037	98362	8577	18325	cm:i:745	NM:i:-1	AS:i:8575	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/7626/19215_23470	4255	29	4225	-	gi|545778205|gb|U00096.3|	100000	94222	98351	1979	4129	cm:i:171	NM:i:-1	AS:i:1977	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/7642/0_6981	6981	317	3537	-	gi|545778205|gb|U00096.3|	100000	119	3443	673	3324	cm:i:61	NM:i:-1	AS:i:673	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/7768/0_20664	20664	14	20477	-	gi|545778205|gb|U00096.3|	100000	78542	97924	10360	19382	cm:i:877	NM:i:-1	AS:i:10351	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/7768/20714_23917	3203	16	3155	-	gi|545778205|gb|U00096.3|	100000	94911	97924	1508	3013	cm:i:130	NM:i:-1	AS:i:1507	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/7998/0_18132	18132	16994	17646	-	gi|545778205|gb|U00096.3|	100000	19797	20534	88	737	cm:i:9	NM:i:-1	AS:i:88	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/8000/0_12334	12334	45	12059	-	gi|545778205|gb|U00096.3|	100000	30187	42417	2981	12230	cm:i:267	NM:i:-1	AS:i:2981	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/8002/0_15731	15731	98	714	-	gi|545778205|gb|U00096.3|	100000	19796	20332	202	536	cm:i:19	NM:i:-1	AS:i:202	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/8133/739_8981	8242	14	8216	-	gi|545778205|gb|U00096.3|	100000	14584	22488	3911	7904	cm:i:342	NM:i:-1	AS:i:3904	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/8156/0_7011	7011	98	6949	-	gi|545778205|gb|U00096.3|	100000	1290	7544	3075	6254	cm:i:259	NM:i:-1	AS:i:3074	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/8426/0_5361	5361	20	1110	-	gi|545778205|gb|U00096.3|	100000	15399	16493	514	1094	cm:i:47	NM:i:-1	AS:i:514	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/8516/0_1610	1610	380	1563	-	gi|545778205|gb|U00096.3|	100000	82445	83579	341	1134	cm:i:27	NM:i:-1	AS:i:341	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/8927/0_18304	18304	12016	13312	-	gi|545778205|gb|U00096.3|	100000	15404	16699	638	1295	cm:i:57	NM:i:-1	AS:i:638	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/8927/18346_35482	17136	11765	13066	-	gi|545778205|gb|U00096.3|	100000	15419	16694	457	1275	cm:i:42	NM:i:-1	AS:i:457	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/8964/0_12026	12026	15	613	-	gi|545778205|gb|U00096.3|	100000	19826	20406	171	580	cm:i:15	NM:i:-1	AS:i:171	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/8976/3577_6456	2879	63	2865	-	gi|545778205|gb|U00096.3|	100000	52591	55442	1013	2851	cm:i:87	NM:i:-1	AS:i:1013	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9104/0_17379	17379	15382	16063	-	gi|545778205|gb|U00096.3|	100000	19836	20537	277	701	cm:i:26	NM:i:-1	AS:i:277	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9158/3746_18052	14306	195	14294	-	gi|545778205|gb|U00096.3|	100000	57994	71219	4002	13225	cm:i:356	NM:i:-1	AS:i:3999	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9164/1351_15175	13824	4196	4971	-	gi|545778205|gb|U00096.3|	100000	19796	20547	450	751	cm:i:37	NM:i:-1	AS:i:449	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9192/0_23164	23164	11064	11820	-	gi|545778205|gb|U00096.3|	100000	19826	20537	381	711	cm:i:31	NM:i:-1	AS:i:381	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9301/0_23092	23092	2045	2809	-	gi|545778205|gb|U00096.3|	100000	19797	20532	468	735	cm:i:41	NM:i:-1	AS:i:468	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9491/18951_24483	5532	5057	5524	-	gi|545778205|gb|U00096.3|	100000	99411	99851	195	440	cm:i:17	NM:i:-1	AS:i:195	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9569/0_3326	3326	293	3288	-	gi|545778205|gb|U00096.3|	100000	40503	43339	1303	2836	cm:i:115	NM:i:-1	AS:i:1300	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9569/3373_6916	3543	16	2926	-	gi|545778205|gb|U00096.3|	100000	40853	43498	1338	2645	cm:i:114	NM:i:-1	AS:i:1337	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9569/6964_10108	3144	32	3101	-	gi|545778205|gb|U00096.3|	100000	40534	43480	1717	2946	cm:i:148	NM:i:-1	AS:i:1717	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9569/10153_13253	3100	53	3056	-	gi|545778205|gb|U00096.3|	100000	40536	43465	1700	2929	cm:i:142	NM:i:-1	AS:i:1700	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9569/13295_16423	3128	40	3094	-	gi|545778205|gb|U00096.3|	100000	40523	43478	1650	2955	cm:i:138	NM:i:-1	AS:i:1648	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9569/16467_19630	3163	17	3159	-	gi|545778205|gb|U00096.3|	100000	40497	43497	1684	3000	cm:i:140	NM:i:-1	AS:i:1684	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9569/19677_22845	3168	17	3164	-	gi|545778205|gb|U00096.3|	100000	40497	43497	1712	3000	cm:i:147	NM:i:-1	AS:i:1710	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9569/22890_25973	3083	44	3059	-	gi|545778205|gb|U00096.3|	100000	40515	43474	1704	2959	cm:i:149	NM:i:-1	AS:i:1700	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9569/26019_29177	3158	156	3090	-	gi|545778205|gb|U00096.3|	100000	40555	43369	1366	2814	cm:i:118	NM:i:-1	AS:i:1365	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9569/29223_32322	3099	17	3053	-	gi|545778205|gb|U00096.3|	100000	40534	43497	1732	2963	cm:i:150	NM:i:-1	AS:i:1729	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9569/32367_33943	1576	27	1562	-	gi|545778205|gb|U00096.3|	100000	42049	43489	774	1440	cm:i:68	NM:i:-1	AS:i:773	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9767/9303_28629	19326	90	19302	-	gi|545778205|gb|U00096.3|	100000	41488	60123	7588	18635	cm:i:657	NM:i:-1	AS:i:7588	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9916/0_15176	15176	341	14031	-	gi|545778205|gb|U00096.3|	100000	3	12918	6541	12915	cm:i:570	NM:i:-1	AS:i:6534	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9916/15222_18611	3389	25	2299	-	gi|545778205|gb|U00096.3|	100000	3	2243	1148	2240	cm:i:101	NM:i:-1	AS:i:1146	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/10301/0_5437	5437	182	5437	-	gi|545778205|gb|U00096.3|	100000	90091	95196	2491	5105	cm:i:215	NM:i:-1	AS:i:2490	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/10307/0_4189	4189	376	4136	-	gi|545778205|gb|U00096.3|	100000	3605	7155	1842	3550	cm:i:155	NM:i:-1	AS:i:1837	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
