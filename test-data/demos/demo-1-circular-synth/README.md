This is a synthetic circular genome, composed of the 100000 bp of the E. Coli genome, in the following manner:
[50000:100000] + [0:50000]

```
echo ">synth-circ-1 Ecoli-50000:100000-0:50000" > synth-circular.fasta
samtools faidx ../ecoli-small/ecoli-0-100000.fasta "gi|545778205|gb|U00096.3|:50000-100000" | grep -v "^>" | tr -d '\n' >> synth-circular.fasta
echo "" >> synth-circular.fasta
samtools faidx ../ecoli-small/ecoli-0-100000.fasta "gi|545778205|gb|U00096.3|:0-50000" | grep -v "^>" | tr -d '\n' >> synth-circular.fasta
echo "" >> synth-circular.fasta
```

Constructing synthetic, exact reads for testing:
```
samtools faidx ../ecoli-small/ecoli-0-100000.fasta "gi|545778205|gb|U00096.3|:40000-59999" > read-1-exact-match-circular.fasta
samtools faidx ../ecoli-small/ecoli-0-100000.fasta "gi|545778205|gb|U00096.3|:60000-79999" > read-2-exact-match-linear.fasta
```

