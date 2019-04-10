#! /bin/bash

# set -vex

if [ "$#" -lt 4 ]; then
  echo "Wrong number of arguments."
  echo ""
  echo "Usage:"
  echo "  mapeval <reads.fasta> <truth.paf> <query.paf> <out_prefix> [run_from_stage]"
  echo ""
  exit
fi

reads=$1
set1=$2
set2=$3
prefix=$4
# This parameter enables to skip some stages, e.g. parsing the read headers.
run_from_stage=1

if [ "$#" -ge 5 ]; then
    run_from_stage=$5
fi

# Set1 is viewed as the "Truth" set. Everything mapped there is considered the full set of
# the input, and everything which wasn't mapped in Set1 is ignored.
# If Set2 contains reads which weren't mapped in Set1, those are discarded.
# Out of the other ones True Positives (TP) are the ones mapped to the same
# location by both Set1 and Set2.
# The False Negatives (FN) are the ones mapped in Set1, but not in Set2.

function parse_reads_headers {
    python -c "import sys; sys.stderr.write('[0/3 - Parsing reads headers]\n')"
    grep ">" ${reads} | tr -d ">" > ${prefix}.input.headers
}

function parse_inputs {
    # Encode qname, rstrand and rname in the ChrName of the header.
    # That way, it gets evaluated as a tuple.
    python -c "import sys; sys.stderr.write('[1/3 - Parsing input PAF files]\n')"
    awk '{ print $1"!"$5"!"$6"\t"$8"\t"$9 }' ${set1} > ${prefix}.set1.bed
    awk '{ print $1"!"$5"!"$6"\t"$8"\t"$9 }' ${set2} > ${prefix}.set2.bed
    awk '{ print $1 }' ${set1} | sort | uniq > ${prefix}.set1.headers
    awk '{ print $1 }' ${set2} | sort | uniq > ${prefix}.set2.headers
    cat ${prefix}.set1.headers > ${prefix}.both.headers
    cat ${prefix}.set2.headers >> ${prefix}.both.headers
    cat ${prefix}.both.headers | sort | uniq -c | sort -k 1,1 -n > ${prefix}.both.headers.count
}

function eval_sets {
    python -c "import sys; sys.stderr.write('[2/3 - Evaluating mappings]\n')"
    bedtools intersect -a ${prefix}.set1.bed -b ${prefix}.set2.bed -f 0.10 -F 0.10  > ${prefix}.sets.intersect.bed
}

function calc_stats {
    python -c "import sys; sys.stderr.write('[3/3 - Calculating stats]\n')"
    set1_count=$(cat ${prefix}.set1.headers | wc -l)
    set2_count=$(cat ${prefix}.set2.headers | wc -l)
    # input_count=$(cat ${prefix}.input.headers | wc -l)
    mappings1_count=$(cat ${prefix}.set1.bed | wc -l)
    mappings2_count=$(cat ${prefix}.set2.bed | wc -l)
    # Split the ChrName column, to deduplicate the qnames.
    tp_count=$(awk '{ split($1, a, "!"); print a[1]}' ${prefix}.sets.intersect.bed | sort | uniq -c | wc -l)

    # How many reads were mapped in both Set1 and Set2
    set_intersect_count=$(awk '$1 != 1 {print $1}' ${prefix}.both.headers.count | wc -l)

    python -c "import sys; sys.stdout.write('Precision Recall TP FP FN in1 in2 in_both mappings1 mappings2\n')"
    python -c "import sys; in1 = ${set1_count}; in2 = ${set2_count}; set_intersect = ${set_intersect_count}; mappings1 = ${mappings1_count}; mappings2 = ${mappings2_count}; TP = ${tp_count}; FN = in1 - set_intersect; FP = set_intersect - TP; prec = 100.0 * float(TP) / float(TP + FP); rec = 100.0 * float(TP) / float(TP + FN); sys.stdout.write('%5.2f%% %5.2f%% %d %d %d %d %d %d %d %d\n' % (prec, rec, TP, FP, FN, in1, in2, set_intersect, mappings1, mappings2));"
}

# Actually not needed at the moment.
# if [ ${run_from_stage} -le 1 ]; then
#     parse_reads_headers
# fi
if [ ${run_from_stage} -le 1 ]; then
    parse_inputs
fi
if [ ${run_from_stage} -le 2 ]; then
    eval_sets
fi
if [ ${run_from_stage} -le 3 ]; then
    calc_stats
fi

python -c "import sys; sys.stderr.write('[Done]\n')"

