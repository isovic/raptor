for fn in reads.pile*.headers
do
    xargs samtools faidx ecoli-m64030_190330_071939.Q20.30x.fasta < ${fn} > ${fn%.headers}.fasta
done
