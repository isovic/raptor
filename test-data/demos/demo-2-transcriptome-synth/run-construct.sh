#! /bin/bash

function get_exons_pathological_1 {
	# The third exon is hard to map.

	read1_fasta=read-1-spliced-linear.fasta
	echo ">spliced_read_1" > ${read1_fasta}

	echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/8133/739_8981	7742	8242" > temp.bed
	bedtools getfasta -fi ../../ecoli-small/reads.6x.fwd.fasta -bed temp.bed | grep -v "^>" | tr -d '\n' >> ${read1_fasta}
	echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9569/16467_19630	0	700" > temp.bed
	bedtools getfasta -fi ../../ecoli-small/reads.6x.fwd.fasta -bed temp.bed | grep -v "^>" | tr -d '\n' >> ${read1_fasta}
	echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9767/9303_28629	0	400" > temp.bed
	bedtools getfasta -fi ../../ecoli-small/reads.6x.fwd.fasta -bed temp.bed | grep -v "^>" | tr -d '\n' >> ${read1_fasta}
	echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/8516/0_1610	610	1610" > temp.bed
	bedtools getfasta -fi ../../ecoli-small/reads.6x.fwd.fasta -bed temp.bed | grep -v "^>" | tr -d '\n' >> ${read1_fasta}
}

function get_exons_linear {

	read1_fasta=read-1-spliced-linear.fasta
	echo ">spliced_read_1" > ${read1_fasta}

	echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/8133/739_8981	7742	8242" > temp.bed
	bedtools getfasta -fi ../../ecoli-small/reads.6x.fwd.fasta -bed temp.bed | grep -v "^>" | tr -d '\n' >> ${read1_fasta}

	echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9569/16467_19630	0	700" > temp.bed
	bedtools getfasta -fi ../../ecoli-small/reads.6x.fwd.fasta -bed temp.bed | grep -v "^>" | tr -d '\n' >> ${read1_fasta}

	echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9767/9303_28629	19000	19325" > temp.bed
	bedtools getfasta -fi ../../ecoli-small/reads.6x.fwd.fasta -bed temp.bed | grep -v "^>" | tr -d '\n' >> ${read1_fasta}

	echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/8516/0_1610	610	1610" > temp.bed
	bedtools getfasta -fi ../../ecoli-small/reads.6x.fwd.fasta -bed temp.bed | grep -v "^>" | tr -d '\n' >> ${read1_fasta}

	# echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/8133/739_8981	7742	8242" > exons.bed
	# echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9569/16467_19630	0	700" >> exons.bed
	# echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9767/9303_28629	0	400" >> exons.bed
	# echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/8516/0_1610	610	1610" >> exons.bed

	# m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/8133/739_8981    8242    25      8242    +       gi|545778205|gb|U00096.3|       100000  14584   22502   8217    7918    60      cm:i:7627
	# m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9569/16467_19630 3163    0       3163    +       gi|545778205|gb|U00096.3|       100000  40493   43514   3163    3021    60      cm:i:2912
	# m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9767/9303_28629  19326   0       19325   +       gi|545778205|gb|U00096.3|       100000  41467   60192   19325   18725   60      cm:i:17746
	# m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/8516/0_1610      1610    45      1610    +       gi|545778205|gb|U00096.3|       100000  82442   83896   1565    1454    60      cm:i:1337
}

function get_exons_circular {
	# Simulate a circRNA by rotating a couple of exons.

	read1_fasta=read-2-spliced-circular.fasta
	echo ">spliced_read_2" > ${read1_fasta}

	echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9767/9303_28629	19000	19325" > temp.bed
	bedtools getfasta -fi ../../ecoli-small/reads.6x.fwd.fasta -bed temp.bed | grep -v "^>" | tr -d '\n' >> ${read1_fasta}

	echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/8516/0_1610	610	1610" > temp.bed
	bedtools getfasta -fi ../../ecoli-small/reads.6x.fwd.fasta -bed temp.bed | grep -v "^>" | tr -d '\n' >> ${read1_fasta}

	echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/8133/739_8981	7742	8242" > temp.bed
	bedtools getfasta -fi ../../ecoli-small/reads.6x.fwd.fasta -bed temp.bed | grep -v "^>" | tr -d '\n' >> ${read1_fasta}

	echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/9569/16467_19630	0	700" > temp.bed
	bedtools getfasta -fi ../../ecoli-small/reads.6x.fwd.fasta -bed temp.bed | grep -v "^>" | tr -d '\n' >> ${read1_fasta}
}

get_exons_linear
get_exons_circular
