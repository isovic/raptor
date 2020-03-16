##############
### Test-1 ###
##############
Run a small test sample with missing adapter. One piece should be reported as a secondary and one as primary.
  $ ${BIN_DIR}/raptor -r ${PROJECT_DIR}/test-data/secondary-suppl/test-1-adapter/ref.fasta -q ${PROJECT_DIR}/test-data/secondary-suppl/test-1-adapter/reads.fasta -v 0 | sed 's/[[:space:]]pi:.*ps:i:[0-9]//g'
  m54081_181222_131903/18219648/54727_74589	19862	18	9745	-	ctg.s1.000000F:53800-63300	9501	143	9414	4914	9271	3	cm:i:415	NM:i:-1	AS:i:4914	fg:i:16	cg:Z:*
  m54081_181222_131903/18219648/54727_74589	19862	10050	19845	+	ctg.s1.000000F:53800-63300	9501	82	9416	4841	9334	3	cm:i:408	NM:i:-1	AS:i:4841	fg:i:256	cg:Z:*

The same test example as above, but testing filtering on mapped length. The primary piece has mapped length 9727 and secondary is larger with a span of 9795. If the primary piece is filtered out by length, the secondary shouldn't be reported either.
  $ ${BIN_DIR}/raptor -r ${PROJECT_DIR}/test-data/secondary-suppl/test-1-adapter/ref.fasta -q ${PROJECT_DIR}/test-data/secondary-suppl/test-1-adapter/reads.fasta -v 0 --min-map-len 9728 | sed 's/[[:space:]]pi:.*ps:i:[0-9]//g'

Test the bestn filter. Setting "--bestn 1" should report only the primary alignment.
  $ ${BIN_DIR}/raptor -r ${PROJECT_DIR}/test-data/secondary-suppl/test-1-adapter/ref.fasta -q ${PROJECT_DIR}/test-data/secondary-suppl/test-1-adapter/reads.fasta -v 0 --bestn 1 --bestn-threshold 1.0 | sed 's/[[:space:]]pi:.*ps:i:[0-9]//g'
  m54081_181222_131903/18219648/54727_74589	19862	18	9745	-	ctg.s1.000000F:53800-63300	9501	143	9414	4914	9271	3	cm:i:415	NM:i:-1	AS:i:4914	fg:i:16	cg:Z:*

Test the bestn filter. Setting "--bestn 10" should report at most 9 secondary alignments and 1 primary alignment. In this case, there is 1 secondary alignment.
  $ ${BIN_DIR}/raptor -r ${PROJECT_DIR}/test-data/secondary-suppl/test-1-adapter/ref.fasta -q ${PROJECT_DIR}/test-data/secondary-suppl/test-1-adapter/reads.fasta -v 0 --bestn 10 --bestn-threshold 1.0 | sed 's/[[:space:]]pi:.*ps:i:[0-9]//g'
  m54081_181222_131903/18219648/54727_74589	19862	18	9745	-	ctg.s1.000000F:53800-63300	9501	143	9414	4914	9271	3	cm:i:415	NM:i:-1	AS:i:4914	fg:i:16	cg:Z:*
  m54081_181222_131903/18219648/54727_74589	19862	10050	19845	+	ctg.s1.000000F:53800-63300	9501	82	9416	4841	9334	3	cm:i:408	NM:i:-1	AS:i:4841	fg:i:256	cg:Z:*

Test the bestn filter. Setting "--bestn 0" should report all secondary alignments and 1 primary alignment. In this case, there is 1 secondary alignment.
  $ ${BIN_DIR}/raptor -r ${PROJECT_DIR}/test-data/secondary-suppl/test-1-adapter/ref.fasta -q ${PROJECT_DIR}/test-data/secondary-suppl/test-1-adapter/reads.fasta -v 0 --bestn 0 --bestn-threshold 1.0 | sed 's/[[:space:]]pi:.*ps:i:[0-9]//g'
  m54081_181222_131903/18219648/54727_74589	19862	18	9745	-	ctg.s1.000000F:53800-63300	9501	143	9414	4914	9271	3	cm:i:415	NM:i:-1	AS:i:4914	fg:i:16	cg:Z:*
  m54081_181222_131903/18219648/54727_74589	19862	10050	19845	+	ctg.s1.000000F:53800-63300	9501	82	9416	4841	9334	3	cm:i:408	NM:i:-1	AS:i:4841	fg:i:256	cg:Z:*

##############
### Test-2 ###
##############
Test a combination of secondary and supplementary mappings on real data.
  $ ${BIN_DIR}/raptor -r ${PROJECT_DIR}/test-data/secondary-suppl/test-2-diff-contig/ref.fasta -q ${PROJECT_DIR}/test-data/secondary-suppl/test-2-diff-contig/reads.fasta -v 0 --align | sed "s/[[:space:]]cg:Z:.*//g"
  m54047_190612_045721/5177668/15116_33236	18120	4415	18120	+	ctg.s2.000009F	37758	0	13162	13705	13162	3	cm:i:12364	NM:i:1796	AS:i:18256	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0
  m54047_190612_045721/5177668/15116_33236	18120	0	727	-	ctg.s2.000001F	59739	8837	9555	727	718	3	cm:i:629	NM:i:131	AS:i:766	fg:i:2064	pi:i:3	pj:i:0	pn:i:1	ps:i:0
  m54047_190612_045721/5177668/15116_33236	18120	0	13750	-	ctg.s2.000009F	37758	31	13162	13750	13131	2	cm:i:12351	NM:i:1833	AS:i:18100	fg:i:272	pi:i:1	pj:i:0	pn:i:1	ps:i:0
  m54047_190612_045721/5177668/33312_51083	17771	4358	17770	+	ctg.s2.000009F	37758	0	13159	13412	13159	2	cm:i:12476	NM:i:1363	AS:i:19892	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0
  m54047_190612_045721/5177668/33312_51083	17771	0	780	-	ctg.s2.000001F	59739	8838	9555	780	717	1	cm:i:642	NM:i:158	AS:i:738	fg:i:2064	pi:i:3	pj:i:0	pn:i:1	ps:i:0
  m54047_190612_045721/5177668/33312_51083	17771	0	13389	-	ctg.s2.000009F	37758	11	13162	13389	13151	1	cm:i:12458	NM:i:1364	AS:i:19780	fg:i:272	pi:i:1	pj:i:0	pn:i:1	ps:i:0
  m54047_190612_045721/5177668/33312_51083	17771	29	794	+	ctg.s2.000001F	59739	16353	17060	765	707	2	cm:i:623	NM:i:164	AS:i:674	fg:i:256	pi:i:2	pj:i:0	pn:i:1	ps:i:0
  m54047_190612_045721/5177668/33312_51083	17771	16963	17770	+	ctg.s2.000001F	59739	8842	9552	807	710	2	cm:i:635	NM:i:197	AS:i:662	fg:i:256	pi:i:4	pj:i:0	pn:i:1	ps:i:0
  m54047_190612_045721/5177668/51159_68911	17752	4332	17752	+	ctg.s2.000009F	37758	9	13162	13420	13153	2	cm:i:12505	NM:i:1328	AS:i:20076	fg:i:0	pi:i:0	pj:i:0	pn:i:1	ps:i:0
  m54047_190612_045721/5177668/51159_68911	17752	1	787	+	ctg.s2.000001F	59739	16324	17099	786	775	1	cm:i:681	NM:i:136	AS:i:834	fg:i:2048	pi:i:3	pj:i:0	pn:i:1	ps:i:0
  m54047_190612_045721/5177668/51159_68911	17752	1	13302	-	ctg.s2.000009F	37758	123	13162	13301	13039	1	cm:i:12408	NM:i:1299	AS:i:19988	fg:i:272	pi:i:1	pj:i:0	pn:i:1	ps:i:0
  m54047_190612_045721/5177668/51159_68911	17752	1	687	-	ctg.s2.000001F	59739	8884	9555	686	671	2	cm:i:609	NM:i:100	AS:i:832	fg:i:272	pi:i:4	pj:i:0	pn:i:1	ps:i:0
  m54047_190612_045721/5177668/51159_68911	17752	16911	17737	-	ctg.s2.000001F	59739	16337	17133	826	796	2	cm:i:692	NM:i:155	AS:i:770	fg:i:272	pi:i:2	pj:i:0	pn:i:1	ps:i:0

##############
### Test-3 ###
##############
Simple synthetic test for supplementary mapping. There are three regions of the E. Coli, permuted. One should be primary and two supplementary.
  $ ${BIN_DIR}/raptor -r ${PROJECT_DIR}/test-data/secondary-suppl/test-3-supp_3_regions/ecoli-0-100000.fasta -q ${PROJECT_DIR}/test-data/secondary-suppl/test-3-supp_3_regions/read-1.fasta -v 0
  gi|545778205|gb|U00096.3|:1-5000	15002	5000	9984	+	gi|545778205|gb|U00096.3|	100000	69999	74983	4998	4984	60	cm:i:362	NM:i:-1	AS:i:4998	fg:i:0	pi:i:2	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  gi|545778205|gb|U00096.3|:1-5000	15002	10001	14984	+	gi|545778205|gb|U00096.3|	100000	29999	34982	4997	4983	60	cm:i:363	NM:i:-1	AS:i:4997	fg:i:2048	pi:i:1	pj:i:0	pn:i:1	ps:i:0	cg:Z:*
  gi|545778205|gb|U00096.3|:1-5000	15002	0	4982	+	gi|545778205|gb|U00096.3|	100000	0	4982	4996	4982	60	cm:i:366	NM:i:-1	AS:i:4996	fg:i:2048	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:*

Same simple synthetic test for supplementary mapping, but with alignment.
  $ ${BIN_DIR}/raptor -r ${PROJECT_DIR}/test-data/secondary-suppl/test-3-supp_3_regions/ecoli-0-100000.fasta -q ${PROJECT_DIR}/test-data/secondary-suppl/test-3-supp_3_regions/read-1.fasta -v 0 --align --aligner ksw2-double
  gi|545778205|gb|U00096.3|:1-5000	15002	10001	15002	+	gi|545778205|gb|U00096.3|	100000	29999	35000	5001	5001	60	cm:i:5001	NM:i:0	AS:i:10002	fg:i:0	pi:i:1	pj:i:0	pn:i:1	ps:i:0	cg:Z:5001=
  gi|545778205|gb|U00096.3|:1-5000	15002	5000	10001	+	gi|545778205|gb|U00096.3|	100000	69999	75000	5001	5001	60	cm:i:5001	NM:i:0	AS:i:10002	fg:i:2048	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:5001=
  gi|545778205|gb|U00096.3|:1-5000	15002	0	5000	+	gi|545778205|gb|U00096.3|	100000	0	5000	5000	5000	60	cm:i:5000	NM:i:0	AS:i:10000	fg:i:2048	pi:i:2	pj:i:0	pn:i:1	ps:i:0	cg:Z:5000=

##############
### Test-4 ###
##############
Synthetic example. Two pieces of a query, one is primary, the other supplementary.
  $ ${BIN_DIR}/raptor -r ${PROJECT_DIR}/test-data/secondary-suppl/test-4-synth-sec_supp_linear/case-1.ref.fasta -q ${PROJECT_DIR}/test-data/secondary-suppl/test-4-synth-sec_supp_linear/case-1.reads.fasta -v 0 --align --aligner ksw2-double
  read-1	10000	0	5000	+	ref	10000	5000	10000	5000	5000	60	cm:i:5000	NM:i:0	AS:i:10000	fg:i:0	pi:i:1	pj:i:0	pn:i:1	ps:i:0	cg:Z:5000=
  read-1	10000	5000	10000	+	ref	10000	0	5000	5000	5000	60	cm:i:5000	NM:i:0	AS:i:10000	fg:i:2048	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:5000=

Synthetic example. Two pieces of a query, one is primary, the other supplementary. The primary portion also has a secondary alignment.
  $ ${BIN_DIR}/raptor -r ${PROJECT_DIR}/test-data/secondary-suppl/test-4-synth-sec_supp_linear/case-2.ref.fasta -q ${PROJECT_DIR}/test-data/secondary-suppl/test-4-synth-sec_supp_linear/case-2.reads.fasta -v 0 --align --aligner ksw2-double
  read-1	10000	0	5000	+	ref	25000	20000	25000	5000	5000	3	cm:i:5000	NM:i:0	AS:i:10000	fg:i:0	pi:i:2	pj:i:0	pn:i:1	ps:i:0	cg:Z:5000=
  read-1	10000	5000	10000	+	ref	25000	0	5000	5000	5000	60	cm:i:5000	NM:i:0	AS:i:10000	fg:i:2048	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:5000=
  read-1	10000	0	5000	+	ref	25000	5000	10000	5000	5000	3	cm:i:5000	NM:i:0	AS:i:10000	fg:i:256	pi:i:1	pj:i:0	pn:i:1	ps:i:0	cg:Z:5000=

Synthetic example. Two pieces of a query, one is primary, the other supplementary. The supplementary also has a secondary alignment.
  $ ${BIN_DIR}/raptor -r ${PROJECT_DIR}/test-data/secondary-suppl/test-4-synth-sec_supp_linear/case-3.ref.fasta -q ${PROJECT_DIR}/test-data/secondary-suppl/test-4-synth-sec_supp_linear/case-3.reads.fasta -v 0 --align --aligner ksw2-double
  read-1	10000	0	5000	+	ref	25000	5000	10000	5000	5000	60	cm:i:5000	NM:i:0	AS:i:10000	fg:i:0	pi:i:2	pj:i:0	pn:i:1	ps:i:0	cg:Z:5000=
  read-1	10000	5000	10000	+	ref	25000	20000	25000	5000	5000	3	cm:i:5000	NM:i:0	AS:i:10000	fg:i:2048	pi:i:1	pj:i:0	pn:i:1	ps:i:0	cg:Z:5000=
  read-1	10000	5000	10000	+	ref	25000	0	5000	5000	5000	3	cm:i:5000	NM:i:0	AS:i:10000	fg:i:256	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:5000=

Synthetic example, graph-based mapping with primary, supplementary and secondary mappings. The graph consists of three nodes. The query has 4 adjacent mapped portions: Part1 maps to the end of "ref-a", Part2 to the beginning of "ref-b" (and covers whole of "ref-b"), Part3 maps to the  beginning of the third node "ref-c", and Part4 is slightly further in the "ref-c" node. Part1 also has a secondary alignment to "ref-c".
  $ ${BIN_DIR}/raptor -r ${PROJECT_DIR}/test-data/secondary-suppl/test-4-synth-sec_supp_linear/case-4.ref.fasta -g ${PROJECT_DIR}/test-data/secondary-suppl/test-4-synth-sec_supp_linear/case-4.graph.gfa2 -q ${PROJECT_DIR}/test-data/secondary-suppl/test-4-synth-sec_supp_linear/case-4.reads.fasta -v 0 --align --aligner ksw2-double
  read-1	20000	0	5000	+	ref-a	10000	5000	10000	5000	5000	3	cm:i:5000	NM:i:0	AS:i:10000	fg:i:0	pi:i:0	pj:i:0	pn:i:3	ps:i:1	cg:Z:5000=
  read-1	20000	5000	10000	+	ref-b	5000	0	5000	5000	5000	60	cm:i:5000	NM:i:0	AS:i:10000	fg:i:2048	pi:i:0	pj:i:1	pn:i:3	ps:i:1	cg:Z:5000=
  read-1	20000	10000	15000	+	ref-c	25000	0	5000	5000	5000	60	cm:i:5000	NM:i:0	AS:i:10000	fg:i:2048	pi:i:0	pj:i:2	pn:i:3	ps:i:1	cg:Z:5000=
  read-1	20000	15000	20000	+	ref-c	25000	10000	15000	5000	5000	60	cm:i:5000	NM:i:0	AS:i:10000	fg:i:2048	pi:i:1	pj:i:0	pn:i:1	ps:i:0	cg:Z:5000=
  read-1	20000	0	5000	+	ref-c	25000	15000	20000	5000	5000	3	cm:i:5000	NM:i:0	AS:i:10000	fg:i:256	pi:i:2	pj:i:0	pn:i:1	ps:i:0	cg:Z:5000=

######################
### Real-1-plasmid ###
######################
Graph alignment of the reference sequences on a small plasmid contig set.
This test is a copy from test-2-circular-aln.t, but is added here for documentation.
It's interesting because there are small segmental duplication with lower identity present in this contig. Those are detected by Raptor, but filtered out with the --sec-ratio 0.8' default.
  $ ${BIN_DIR}/raptor --out-fmt m4 -r ${PROJECT_DIR}/test-data/graph-mapping/real-1-plasmid/asm.plasmids.fa -g ${PROJECT_DIR}/test-data/graph-mapping/real-1-plasmid/contig.s2.gfa2 -q ${PROJECT_DIR}/test-data/graph-mapping/real-1-plasmid/ref.plasmids.fa -v 0 --align
  gi|386611788|ref|NC_017637.1| ctg.s2.000000F -25018 99.7935 0 0 12593 102536 1 0 12571 102273
  gi|386611788|ref|NC_017637.1| ctg.s2.000000F -177984 99.6731 0 12593 102536 102536 1 12571 102273 102273
  gi|386607294|ref|NC_017636.1| 1 -7798 99.2974 0 0 3985 5360 0 1372 5332 5332
  gi|386607294|ref|NC_017636.1| 1 -2716 99.6364 0 3985 5360 5360 0 0 1372 5332
# These are removed because of the secondary-to-primary score ratio (default of 0.8):
#  gi|386611788|ref|NC_017637.1| ctg.s2.000000F -644 86.8231 0 6484 7038 102536 1 11903 12439 102273
#  gi|386611788|ref|NC_017637.1| ctg.s2.000000F -628 89.8488 0 132 595 102536 1 5648 6101 102273
#  gi|386611788|ref|NC_017637.1| ctg.s2.000000F -518 89.4188 0 5403 5781 102536 1 13738 14099 102273
#  gi|386611788|ref|NC_017637.1| ctg.s2.000000F -486 95.7597 0 101080 101362 102536 1 6806 7089 102273

Graph alignment of the reference sequences on a small plasmid contig set.
The same as above, but the secondary-to-primary score ratio was set to 0 with '--sec-ratio 0.0'. This allows for detection of lower-quality segmental duplications.
  $ ${BIN_DIR}/raptor --out-fmt m4 -r ${PROJECT_DIR}/test-data/graph-mapping/real-1-plasmid/asm.plasmids.fa -g ${PROJECT_DIR}/test-data/graph-mapping/real-1-plasmid/contig.s2.gfa2 -q ${PROJECT_DIR}/test-data/graph-mapping/real-1-plasmid/ref.plasmids.fa -v 0 --align --sec-ratio 0.0
  gi|386611788|ref|NC_017637.1| ctg.s2.000000F -25018 99.7935 0 0 12593 102536 1 0 12571 102273
  gi|386611788|ref|NC_017637.1| ctg.s2.000000F -177984 99.6731 0 12593 102536 102536 1 12571 102273 102273
  gi|386611788|ref|NC_017637.1| ctg.s2.000000F -644 86.8231 0 6484 7038 102536 1 11903 12439 102273
  gi|386611788|ref|NC_017637.1| ctg.s2.000000F -628 89.8488 0 132 595 102536 1 5648 6101 102273
  gi|386611788|ref|NC_017637.1| ctg.s2.000000F -518 89.4180 0 5403 5781 102536 1 13738 14099 102273
  gi|386611788|ref|NC_017637.1| ctg.s2.000000F -486 95.7597 0 101080 101362 102536 1 6806 7089 102273
  gi|386607294|ref|NC_017636.1| 1 -7798 99.2974 0 0 3985 5360 0 1372 5332 5332
  gi|386607294|ref|NC_017636.1| 1 -2716 99.6364 0 3985 5360 5360 0 0 1372 5332



##############################
### Mapping quality tests. ###
##############################
A read which maps only partially to the reference. The mapped portion is very small, and even though it's unique, the mapping quality will be >3 and <60.
  $ ${BIN_DIR}/raptor --align -n 100 -r ${PROJECT_DIR}/test-data/ecoli-small/ecoli-0-100000.fasta -q ${PROJECT_DIR}/test-data/ecoli-small/reads.6x.rev.fasta -v 0 -n 1 -s 1 | sed -E 's/AS:i:[0-9]+[[:space:]]//g'
  m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3414/0_11983	11983	164	1053	-	gi|545778205|gb|U00096.3|	100000	19788	20565	889	777	4	cm:i:722	NM:i:204	fg:i:16	pi:i:0	pj:i:0	pn:i:1	ps:i:0	cg:Z:1D1X2=1I1X1=1X8=1I31=1I6=2D10=1I5=1I2=1D3=1I3=1D22=1I9=1D1=1X14=1I2=1I15=1I3=2I3=1D12=2I4=1D18=1D4=1D8=1I8=1D8=1D7=1D13=1D1=1X13=1I17=1X1=2D9=1D16=1D9=2D5=1I4=1D1=1D5=1X7=3I1=1I20=2I7=1I12=1X2=1D1X4=1D4=1D5=1D2=1D10=1I3=1D8=1D11=1I2=1I24=1I8=1X2=1X2=1X15=1D5=1D4=1I4=1I5=1X10=1I10=1I15=1D26=2I11=1D15=1D4=3I14=1I20=1D5=1D9=1I19=1I1=3I2=1I1=1X1=1I1=1I1=2I1=4I2=1I2=1I2=2I1=4I3=7I2=6I3=3I2=4I1=1I3=11I4=1I2=14I3=8I2=2I2=3I1=10I2=7I2=3I4=1I11=1I4=1I17=1I4=1I1=2I1=2I1=1X2=3X1=1D1=


##############################
### Filtering by identity. ###
##############################
The following sequence in asm.plasmids.fa aligns circularly to the reference. The longer piece (primary alignment) has slightly lower identity.
If the `--min-idt` filter is applied so that the primary alignment gets filtered out, then what was supplementary alignment now becomes the primary.
  $ ${BIN_DIR}/raptor --out-fmt m4 -r ${PROJECT_DIR}/test-data/graph-mapping/real-1-plasmid/ref.plasmids.fa -q ${PROJECT_DIR}/test-data/graph-mapping/real-1-plasmid/asm.plasmids.fa -v 0 --align --out-fmt paf -n 1 -s 0 --min-idt 99.68
  ctg.s2.000000F	102273	0	12571	-	gi|386611788|ref|NC_017637.1|	102536	0	12593	12571	12593	4	cm:i:12567	NM:i:30	AS:i:25018	fg:i:16	pi:i:1	pj:i:0	pn:i:1	ps:i:0	cg:Z:248=1D348=1D194=1D186=1D566=1D1518=1D501=1D217=1D1704=1D61=1D220=1D208=1D160=1D90=1D280=1D1898=1D111=1D442=1D455=1D178=2I84=1I233=1D343=1D573=1D93=1D199=1D158=2D910=1I389=
