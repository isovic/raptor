#! /bin/bash

function make_ref {
	echo ">synth-circ-1 Ecoli-50000:100000-0:50000" > synth-circular.fasta
	echo "gi|545778205|gb|U00096.3|	50000	100000" > temp.bed
	bedtools getfasta -fi ../ecoli-small/ecoli-0-100000.fasta -bed temp.bed | grep -v "^>" | tr -d '\n' >> synth-circular.fasta
	echo "" >> synth-circular.fasta
        echo "gi|545778205|gb|U00096.3|	0	50000" > temp.bed
        bedtools getfasta -fi ../ecoli-small/ecoli-0-100000.fasta -bed temp.bed | grep -v "^>" | tr -d '\n' >> synth-circular.fasta
	echo "" >> synth-circular.fasta
}

function make_reads {
        echo "gi|545778205|gb|U00096.3|	40000	60000" > temp.bed
        bedtools getfasta -fi ../ecoli-small/ecoli-0-100000.fasta -bed temp.bed > read-1-exact-match-circular.fasta

        echo "gi|545778205|gb|U00096.3|	60000	80000" > temp.bed
        bedtools getfasta -fi ../ecoli-small/ecoli-0-100000.fasta -bed temp.bed > read-2-exact-match-linear.fasta
}

make_ref
make_reads

