#! /bin/bash

function get_tiles {
  samtools faidx ecoli-0-100000.fasta "gi|545778205|gb|U00096.3|:30001-35000" > part.ref.1.30k-35k.fasta
  samtools faidx ecoli-0-100000.fasta "gi|545778205|gb|U00096.3|:35001-40000" > part.ref.2.35k-40k.fasta
  samtools faidx ecoli-0-100000.fasta "gi|545778205|gb|U00096.3|:40001-45000" > part.ref.3.40k-45k.fasta
  samtools faidx ecoli-0-100000.fasta "gi|545778205|gb|U00096.3|:45001-50000" > part.ref.4.45k-50k.fasta
  samtools faidx ecoli-0-100000.fasta "gi|545778205|gb|U00096.3|:50001-55000" > part.ref.5.50k-55k.fasta
  samtools faidx ecoli-0-100000.fasta "gi|545778205|gb|U00096.3|:55001-60000" > part.ref.6.55k-60k.fasta
  samtools faidx ecoli-0-100000.fasta "gi|545778205|gb|U00096.3|:60001-65000" > part.ref.7.60k-65k.fasta
  samtools faidx ecoli-0-100000.fasta "gi|545778205|gb|U00096.3|:65001-70000" > part.ref.8.65k-70k.fasta
  samtools faidx ecoli-0-100000.fasta "gi|545778205|gb|U00096.3|:70001-75000" > part.ref.9.70k-75k.fasta
  samtools faidx ecoli-0-100000.fasta "gi|545778205|gb|U00096.3|:75001-80000" > part.ref.10.75k-80k.fasta
}

TILE1=part.ref.1.30k-35k.fasta
TILE2=part.ref.2.35k-40k.fasta
TILE3=part.ref.3.40k-45k.fasta
TILE4=part.ref.4.45k-50k.fasta
TILE5=part.ref.5.50k-55k.fasta
TILE6=part.ref.6.55k-60k.fasta
TILE7=part.ref.7.60k-65k.fasta
TILE8=part.ref.8.65k-70k.fasta
TILE9=part.ref.9.70k-75k.fasta
TILE10=part.ref.10.75k-80k.fasta

function combine {
    local local_out_fn=$1
    shift
    local local_header=$1
    shift
    local local_tiles=$@
    echo ">${local_header}" >> ${local_out_fn}
    for tile in ${local_tiles}
    do
        tail -n +2 ${tile} >> ${local_out_fn}
    done
}

function test_case_1 {
    # Primary and secondary alignments.
    local local_id=$1
    local local_out_ref="case-${local_id}.ref.fasta"
    local local_out_reads="case-${local_id}.reads.fasta"
    rm -f ${local_out_ref}
    rm -f ${local_out_reads}
    combine ${local_out_ref} "ref" ${TILE2} ${TILE1}
    combine ${local_out_reads} "read-1" ${TILE1} ${TILE2}
}

function test_case_2 {
    # The primary portion also has a secondary alignment.
    local local_id=$1
    local local_out_ref="case-${local_id}.ref.fasta"
    local local_out_reads="case-${local_id}.reads.fasta"
    rm -f ${local_out_ref}
    rm -f ${local_out_reads}
    combine ${local_out_ref} "ref" ${TILE2} ${TILE1} ${TILE3} ${TILE4} ${TILE1}
    combine ${local_out_reads} "read-1" ${TILE1} ${TILE2}
}

function test_case_3 {
    # The supplementary portion has a secondary alignment.
    local local_id=$1
    local local_out_ref="case-${local_id}.ref.fasta"
    local local_out_reads="case-${local_id}.reads.fasta"
    rm -f ${local_out_ref}
    rm -f ${local_out_reads}
    combine ${local_out_ref} "ref" ${TILE2} ${TILE1} ${TILE3} ${TILE4} ${TILE2}
    combine ${local_out_reads} "read-1" ${TILE1} ${TILE2}
}

function test_case_4 {
    # Graph-based test dataset for a combination of secondary/supplemenetary mappings.
    # There are three nodes in a linear graph, the first node has the first match,
    # second node the second match, and the third node the third match plus a repeat
    # of the first match.

    local local_id=$1
    local local_out_ref="case-${local_id}.ref.fasta"
    local local_out_reads="case-${local_id}.reads.fasta"
    local local_out_graph="case-${local_id}.graph.gfa2"
    rm -f ${local_out_ref}
    rm -f ${local_out_reads}

    # Make the reference sequences.
    combine "case-${local_id}.ref-a.fasta" "ref-a" ${TILE4} ${TILE1}
    combine "case-${local_id}.ref-b.fasta" "ref-b" ${TILE2}
    combine "case-${local_id}.ref-c.fasta" "ref-c" ${TILE3} ${TILE5} ${TILE6} ${TILE1} ${TILE7}
    # combine "case-${local_id}.ref-b.fasta" "ref-b" ${TILE2} ${TILE5}
    # combine "case-${local_id}.ref-c.fasta" "ref-c" ${TILE6} ${TILE3} ${TILE7} ${TILE1} ${TILE8}
    cat "case-${local_id}.ref-a.fasta" "case-${local_id}.ref-b.fasta" "case-${local_id}.ref-c.fasta" > ${local_out_ref}
    rm "case-${local_id}.ref-a.fasta" "case-${local_id}.ref-b.fasta" "case-${local_id}.ref-c.fasta"

    # Make the graph.
    echo "H	VN:Z:2.0" > ${local_out_graph}
    echo "S	ref-a	5000	*" >> ${local_out_graph}
    echo "S	ref-b	5000	*" >> ${local_out_graph}
    echo "S	ref-c	25000	*" >> ${local_out_graph}
    echo "E	edge-ab	ref-a+	ref-b+	10000$	10000$	0	0	*" >> ${local_out_graph}
    echo "E	edge-bc	ref-b+	ref-c+	5000$	5000$	0	0	*" >> ${local_out_graph}

    # Make the query sequences.
    combine ${local_out_reads} "read-1" ${TILE1} ${TILE2} ${TILE3} ${TILE6}
    # combine ${local_out_reads} "read-1" ${TILE1} ${TILE2} ${TILE3}
}

# test_case_1 1
# test_case_2 2
# test_case_3 3
test_case_4 4

