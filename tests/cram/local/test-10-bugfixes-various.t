Removing the identity filter exposed the problem with calculating the Alignment Score in case of aligner_edlib. The AS was simply the number of matches, which turned out to be wrong. If only number of matches is used, then indels are ignored, and lower quality alignments can prevail
The path ID information is removed using sed, so that internal labeling does not break the tests.
  $ ${BIN_DIR}/raptor -r ${PROJECT_DIR}/test-data/bugfixes/min_idt_issue/ctg.s1.000000F-95000-115000.fasta -g ${PROJECT_DIR}/test-data/bugfixes/min_idt_issue/all_contig.gfa2 -q ${PROJECT_DIR}/test-data/bugfixes/min_idt_issue/problematic_read.fasta -v 0 --align --out-fmt sam --bestn 0 --bestn-threshold 1.0 --min-idt 0.0 | sed 's/[[:space:]]pi:.*ps:i:[0-9]//g'
  @HD	VN:1.5
  @SQ	SN:ctg.s1.000000F:95000-115000	LN:20001
  m54081_181221_163846/11666180/49398_49554	0	ctg.s1.000000F:95000-115000	13308	46	1=1I70=1X54=2X1I1=1I1=1I1=3X18S	*	0	0	TTCATCTGCGCGGGAATGACGATTCAGAAGTTACACGAAACTCAAAAAAAACGAAACCGAACGAACCGGATTTCCGCTTTTACGGGAATGACGGCGCATAAGTTCCCGTGCGGACAGACCTAGATTCGAGAACACAAAAAAAACCAAAAAGGGGGT	*	NM:i:10	AS:i:216	QS:i:0	QE:i:138	QL:i:156	TS:i:13307	TE:i:13441	TL:i:20001
  m54081_181221_163846/11666180/49398_49554	256	ctg.s1.000000F:95000-115000	4343	46	2S2=1X5=1X21=1X1=1I1=1D5=1I1=1X8=1D12=1X15=2X45=5I1=1I2=20S	*	0	0	TTCATCTGCGCGGGAATGACGATTCAGAAGTTACACGAAACTCAAAAAAAACGAAACCGAACGAACCGGATTTCCGCTTTTACGGGAATGACGGCGCATAAGTTCCCGTGCGGACAGACCTAGATTCGAGAACACAAAAAAAACCAAAAAGGGGGT	*	NM:i:17	AS:i:178	QS:i:2	QE:i:136	QL:i:156	TS:i:4342	TE:i:4470	TL:i:20001
  m54081_181221_163846/11666180/49398_49554	256	ctg.s1.000000F:95000-115000	5151	46	1=1X2=1X27=1X1=1X6=1I1=1X6=1I1=1D12=1X15=2X14=1X25=1X4=2X1I1=1X1=2X1=1X1=2X1=2X13S	*	0	0	TTCATCTGCGCGGGAATGACGATTCAGAAGTTACACGAAACTCAAAAAAAACGAAACCGAACGAACCGGATTTCCGCTTTTACGGGAATGACGGCGCATAAGTTCCCGTGCGGACAGACCTAGATTCGAGAACACAAAAAAAACCAAAAAGGGGGT	*	NM:i:24	AS:i:144	QS:i:0	QE:i:143	QL:i:156	TS:i:5150	TE:i:5291	TL:i:20001

The same as above, but outputting only the primary alignment.
  $ ${BIN_DIR}/raptor -r ${PROJECT_DIR}/test-data/bugfixes/min_idt_issue/ctg.s1.000000F-95000-115000.fasta -g ${PROJECT_DIR}/test-data/bugfixes/min_idt_issue/all_contig.gfa2 -q ${PROJECT_DIR}/test-data/bugfixes/min_idt_issue/problematic_read.fasta -v 0 --align --out-fmt sam | sed 's/[[:space:]]pi:.*ps:i:[0-9]//g'
  @HD	VN:1.5
  @SQ	SN:ctg.s1.000000F:95000-115000	LN:20001
  m54081_181221_163846/11666180/49398_49554	0	ctg.s1.000000F:95000-115000	13308	46	1=1I70=1X54=2X1I1=1I1=1I1=3X18S	*	0	0	TTCATCTGCGCGGGAATGACGATTCAGAAGTTACACGAAACTCAAAAAAAACGAAACCGAACGAACCGGATTTCCGCTTTTACGGGAATGACGGCGCATAAGTTCCCGTGCGGACAGACCTAGATTCGAGAACACAAAAAAAACCAAAAAGGGGGT	*	NM:i:10	AS:i:216	QS:i:0	QE:i:138	QL:i:156	TS:i:13307	TE:i:13441	TL:i:20001

Testing multiple read groups in the input BAM file and using SAM output.
  $ ${BIN_DIR}/raptor -r ${PROJECT_DIR}/test-data/bugfixes/multiple_rg/ecoli-0-100000.fasta -q ${PROJECT_DIR}/test-data/bugfixes/multiple_rg/test_multi_rg.bam --out-fmt sam -v 0 | samtools view -H | sort
  @HD	VN:1.5
  @PG	ID:baz2bam	PN:baz2bam	VN:5.0.0.6236	CL:/opt/pacbio/ppa-5.0.0/bin/baz2bam /data/pa/m54110_171115_145004.baz -o /data/pa/m54110_171115_145004 --metadata /data/pa/.m54110_171115_145004.metadata.xml -j 12 -b 12 --progress --silent --minSubLength 50 --minSnr 3.750000 --adapters /data/pa/m54110_171115_145004.adapters.fasta
  @PG	ID:bazFormat	PN:bazformat	VN:1.3.0
  @PG	ID:bazwriter	PN:bazwriter	VN:5.0.0.6236
  @PG	ID:lima	VN:1.0.0 (commit de983e0)
  @RG	ID:6b04dbdd	PL:PACBIO	DS:READTYPE=SUBREAD;Ipd:CodecV1=ip;PulseWidth:CodecV1=pw;BINDINGKIT=100-862-200;SEQUENCINGKIT=101-309-500;BASECALLERVERSION=5.0.0.6236;FRAMERATEHZ=80.000000;BarcodeFile=/pbi/dept/secondary/siv/smrtlink/smrtlink-beta/smrtsuite_166987/install/smrtlink-release_5.1.0.25817/bundles/smrtinub/current/private/pacbio/barcodes/Sequel_RSII_384_barcodes_v1/Sequel_RSII_384_barcodes_v1.barcodeset.xml;BarcodeHash=979b171f27a290d53f63edf9969c14cc;BarcodeCount=384;BarcodeMode=Symmetric;BarcodeQuality=Score	PU:m54110_171116_010407	PM:SEQUEL
  @RG	ID:8268f187	PL:PACBIO	DS:READTYPE=SUBREAD;Ipd:CodecV1=ip;PulseWidth:CodecV1=pw;BINDINGKIT=100-862-200;SEQUENCINGKIT=101-309-500;BASECALLERVERSION=5.0.0.6236;FRAMERATEHZ=80.000000;BarcodeFile=/pbi/dept/secondary/siv/smrtlink/smrtlink-beta/smrtsuite_166987/install/smrtlink-release_5.1.0.25817/bundles/smrtinub/current/private/pacbio/barcodes/Sequel_RSII_384_barcodes_v1/Sequel_RSII_384_barcodes_v1.barcodeset.xml;BarcodeHash=979b171f27a290d53f63edf9969c14cc;BarcodeCount=384;BarcodeMode=Symmetric;BarcodeQuality=Score	PU:m54110_171115_145004	PM:SEQUEL
  @SQ	SN:gi|545778205|gb|U00096.3|	LN:100000

Testing multiple read groups in the input BAM file and using BAM output.
  $ ${BIN_DIR}/raptor -r ${PROJECT_DIR}/test-data/bugfixes/multiple_rg/ecoli-0-100000.fasta -q ${PROJECT_DIR}/test-data/bugfixes/multiple_rg/test_multi_rg.bam --out-fmt bam -v 0 | samtools view -H | sort
  @HD	VN:1.5	SO:unknown	pb:3.0.7
  @PG	ID:baz2bam	PN:baz2bam	VN:5.0.0.6236	CL:/opt/pacbio/ppa-5.0.0/bin/baz2bam /data/pa/m54110_171115_145004.baz -o /data/pa/m54110_171115_145004 --metadata /data/pa/.m54110_171115_145004.metadata.xml -j 12 -b 12 --progress --silent --minSubLength 50 --minSnr 3.750000 --adapters /data/pa/m54110_171115_145004.adapters.fasta
  @PG	ID:bazFormat	PN:bazformat	VN:1.3.0
  @PG	ID:bazwriter	PN:bazwriter	VN:5.0.0.6236
  @PG	ID:lima	VN:1.0.0 (commit de983e0)
  @RG	ID:6b04dbdd	PL:PACBIO	DS:READTYPE=SUBREAD;Ipd:CodecV1=ip;PulseWidth:CodecV1=pw;BINDINGKIT=100-862-200;SEQUENCINGKIT=101-309-500;BASECALLERVERSION=5.0.0.6236;FRAMERATEHZ=80.000000;BarcodeFile=/pbi/dept/secondary/siv/smrtlink/smrtlink-beta/smrtsuite_166987/install/smrtlink-release_5.1.0.25817/bundles/smrtinub/current/private/pacbio/barcodes/Sequel_RSII_384_barcodes_v1/Sequel_RSII_384_barcodes_v1.barcodeset.xml;BarcodeHash=979b171f27a290d53f63edf9969c14cc;BarcodeCount=384;BarcodeMode=Symmetric;BarcodeQuality=Score	PU:m54110_171116_010407	PM:SEQUEL
  @RG	ID:8268f187	PL:PACBIO	DS:READTYPE=SUBREAD;Ipd:CodecV1=ip;PulseWidth:CodecV1=pw;BINDINGKIT=100-862-200;SEQUENCINGKIT=101-309-500;BASECALLERVERSION=5.0.0.6236;FRAMERATEHZ=80.000000;BarcodeFile=/pbi/dept/secondary/siv/smrtlink/smrtlink-beta/smrtsuite_166987/install/smrtlink-release_5.1.0.25817/bundles/smrtinub/current/private/pacbio/barcodes/Sequel_RSII_384_barcodes_v1/Sequel_RSII_384_barcodes_v1.barcodeset.xml;BarcodeHash=979b171f27a290d53f63edf9969c14cc;BarcodeCount=384;BarcodeMode=Symmetric;BarcodeQuality=Score	PU:m54110_171115_145004	PM:SEQUEL
  @SQ	SN:gi|545778205|gb|U00096.3|	LN:100000
