minimap2 -c ecoli-0-100000.fasta read-1-12kbp_rotated-fwd.fasta
minimap2 -c ecoli-0-100000.fasta read-2-12kbp_rotated-rev.fasta


raptor --align -r ecoli-0-100000.fasta -d read-1-12kbp_rotated-fwd.fasta -g graph.gfa2 -v 0
# read-1-12kbp_rotated-fwd    100000  0   88000   +   gi|545778205|gb|U00096.3|   100000  12000   100000  88000   88000   60  cm:i:88000  nm:i:0  as:i:88000  pi:i:0  pj:i:0  pn:i:2  ps:i:1  cg:Z:88000=
# read-1-12kbp_rotated-fwd    100000  88000   100000  +   gi|545778205|gb|U00096.3|   100000  0   12000   12000   12000   60  cm:i:12000  nm:i:0  as:i:12000  pi:i:0  pj:i:1  pn:i:2  ps:i:1  cg:Z:12000=

raptor --align -r ecoli-0-100000.fasta -d read-2-12kbp_rotated-rev.fasta -g graph.gfa2 -v 0
# read-2-12kbp_rotated-ref    100000  0   12000   -   gi|545778205|gb|U00096.3|   100000  0   12000   12000   12000   60  cm:i:12000  nm:i:0  as:i:12000  pi:i:0  pj:i:0  pn:i:2  ps:i:1  cg:Z:12000=
# read-2-12kbp_rotated-ref    100000  12000   100000  -   gi|545778205|gb|U00096.3|   100000  12000   100000  88000   88000   60  cm:i:88000  nm:i:0  as:i:88000  pi:i:0  pj:i:1  pn:i:2  ps:i:1  cg:Z:88000=

