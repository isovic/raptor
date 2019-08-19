It looks like there is a bug in filtering the alignments based on min-idt.

Running the defaults results in this:
```
$ raptor -r data/collected_ctg.fasta -g data/all_contig.gfa2 -q data/problematic_read.fasta -v 0 --align --out-fmt sam
@HD     VN:1.5
@SQ     SN:ctg.s1.000000F       LN:2190691
m54081_181221_163846/11666180/49398_49554       0       ctg.s1.000000F  108307  39      1=1I70=1X54=2X1I1=1I1=1I1=3X18S *       0       0       TTCATCTGCGCGGGAATGACGATTCAGAAGTTACACGAAACTCAAAAAAAACGAAACCGAACGAACCGGATTTCCGCTTTTACGGGAATGACGGCGCATAAGTTCCCGTGCGGACAGACCTAGATTCGAGAACACAAAAAAAACCAAAAAGGGGGT    *       NM:i:10 AS:i:128        QS:i:0  QE:i:138        QL:i:156        TS:i:108306     TE:i:108440     TL:i:2190691    pi:i:0  pj:i:0  pn:i:1  ps:i:0
```

Specifying the `--min-idt 0.8` results in a weird alignment:
```
$ raptor -r data/collected_ctg.fasta -g data/all_contig.gfa2 -q data/problematic_read.fasta -v 0 --align --out-fmt sam --min-idt 0.8
@HD     VN:1.5
@SQ     SN:ctg.s1.000000F       LN:2190691
m54081_181221_163846/11666180/49398_49554       256     ctg.s1.000000F  99592   30      2S2=1X5=1X21=1X1=1X37=1X21=23D2=6D1=2D1=1X1=5D1=1D1=3D1=2D1=1D4=10D1=3D1=2D2=2D2=2D1=1D1=2D2=2D2=1D1=3D1=3D3=3D1=2D1=1D1=1D1=2D3=4D1=3D1=2D1=13D2=2D1=6D1=1D2=12D1=1D1=14D13S   *      0
0       TTCATCTGCGCGGGAATGACGATTCAGAAGTTACACGAAACTCAAAAAAAACGAAACCGAACGAACCGGATTTCCGCTTTTACGGGAATGACGGCGCATAAGTTCCCGTGCGGACAGACCTAGATTCGAGAACACAAAAAAAACCAAAAAGGGGGT    *       NM:i:147        AS:i:135        QS:i:2  QE:i:143        QL:i:156        TS:i:99591      TE:i:99873      TL:i:2190691    pi:i:1  pj:i:0  pn:i:1  ps:i:0
```

When I list all of the alignments, this is what I get:
```
$ raptor -r data/collected_ctg.fasta -g data/all_contig.gfa2 -q data/problematic_read.fasta -v 0 --align --out-fmt sam --min-qlen 50 --bestn 0 --bestn-threshold 1.0 --min-idt 0.8
@HD     VN:1.5
@SQ     SN:ctg.s1.000000F       LN:2190691
m54081_181221_163846/11666180/49398_49554       256     ctg.s1.000000F  99592   30      2S2=1X5=1X21=1X1=1X37=1X21=23D2=6D1=2D1=1X1=5D1=1D1=3D1=2D1=1D4=10D1=3D1=2D2=2D2=2D1=1D1=2D2=2D2=1D1=3D1=3D3=3D1=2D1=1D1=1D1=2D3=4D1=3D1=2D1=13D2=2D1=6D1=1D2=12D1=1D1=14D13S   *      0
0       TTCATCTGCGCGGGAATGACGATTCAGAAGTTACACGAAACTCAAAAAAAACGAAACCGAACGAACCGGATTTCCGCTTTTACGGGAATGACGGCGCATAAGTTCCCGTGCGGACAGACCTAGATTCGAGAACACAAAAAAAACCAAAAAGGGGGT    *       NM:i:147        AS:i:135        QS:i:2  QE:i:143        QL:i:156        TS:i:99591      TE:i:99873      TL:i:2190691    pi:i:1  pj:i:0  pn:i:1  ps:i:0
m54081_181221_163846/11666180/49398_49554       0       ctg.s1.000000F  108307  30      1=1I70=1X54=2X1I1=1I1=1I1=3X18S *       0       0       TTCATCTGCGCGGGAATGACGATTCAGAAGTTACACGAAACTCAAAAAAAACGAAACCGAACGAACCGGATTTCCGCTTTTACGGGAATGACGGCGCATAAGTTCCCGTGCGGACAGACCTAGATTCGAGAACACAAAAAAAACCAAAAAGGGGGT    *       NM:i:10 AS:i:128        QS:i:0  QE:i:138        QL:i:156        TS:i:108306     TE:i:108440     TL:i:2190691    pi:i:0  pj:i:0  pn:i:1  ps:i:0
m54081_181221_163846/11666180/49398_49554       256     ctg.s1.000000F  100150  30      1=1X2=1X27=1X1=1X6=1I1=1X6=1I1=1D12=1X15=2X14=1X25=1X4=2X1I1=1X1=2X1=1X1=2X1=2X13S      *       0       0       TTCATCTGCGCGGGAATGACGATTCAGAAGTTACACGAAACTCAAAAAAAACGAAACCGAACGAACCGGATTTCCGCTTTTACGGGAATGACGGCGCATAAGTTCCCGTGCGGACAGACCTAGATTCGAGAACACAAAAAAAACCAAAAAGGGGGT    *       NM:i:24 AS:i:120        QS:i:0  QE:i:143        QL:i:156        TS:i:100149     TE:i:100290     TL:i:2190691    pi:i:3  pj:i:0  pn:i:1  ps:i:0
m54081_181221_163846/11666180/49398_49554       256     ctg.s1.000000F  99342   30      2S2=1X5=1X21=1X1=1I1=1D5=1I1=1X8=1D12=1X15=2X45=5I1=1I2=20S     *       0       0       TTCATCTGCGCGGGAATGACGATTCAGAAGTTACACGAAACTCAAAAAAAACGAAACCGAACGAACCGGATTTCCGCTTTTACGGGAATGACGGCGCATAAGTTCCCGTGCGGACAGACCTAGATTCGAGAACACAAAAAAAACCAAAAAGGGGGT    *       NM:i:17 AS:i:119        QS:i:2  QE:i:136        QL:i:156        TS:i:99341      TE:i:99469      TL:i:2190691    pi:i:2  pj:i:0  pn:i:1  ps:i:0
```

# Conclusion:
- It turns out that the alignment with the highest AS score is actually the worst one. Why does this happen? How is this possible?
    - I checked the code. The aligner_edlib uses the number of matches as the "score". This is likely wrong, because it doesn't penalize deletions. We should rescore the CIGAR string given the specified parameters.

- The `--min-idt` setting was wrong - the parameter expects the percentage and not the fraction! This allowed alignments with very bad identity to pass through. What should have been set is `--min-idt 80`.