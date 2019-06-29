samtools faidx ecoli-0-100000.fasta "gi|545778205|gb|U00096.3|:1-12000" > temp1.fasta
samtools faidx ecoli-0-100000.fasta "gi|545778205|gb|U00096.3|:12001-100000" > temp2.fasta
echo ">read-1-12kbp_rotated-fwd" > read-1-12kbp_rotated-fwd.fasta
tail -n +2 temp2.fasta >> read-1-12kbp_rotated-fwd.fasta
tail -n +2 temp1.fasta >> read-1-12kbp_rotated-fwd.fasta

~/Documents/git-my/samscripts/src/fastqfilter.py reverse read-1-12kbp_rotated-fwd.fasta > temp3.fasta
echo ">read-2-12kbp_rotated-ref" > read-2-12kbp_rotated-rev.fasta
tail -n +2 temp3.fasta >> read-2-12kbp_rotated-rev.fasta

rm temp1.fasta temp2.fasta temp3.fasta

