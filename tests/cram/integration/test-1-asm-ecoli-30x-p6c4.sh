# SCRIPT_PATH="$( cd "$(dirname "$0")" ; pwd -P )"
# PROJECT_DIR=${SCRIPT_PATH}/../../../

reads=${PROJECT_DIR}/raptor-test-data/overlap/ecoli-30x-p6c4/reads.fasta.gz
paf=out.test.paf
gfa=out.test.gfa

# Overlap.
${BIN_DIR}/raptor -r ${reads} -d ${reads} --overlapper -B 1000 --min-map-len 1000 --bestn 0 --bestn-threshold 1.0 -t 4 -o out.test.paf -v 0

# Assemble.
${PROJECT_DIR}/tools/miniasm/miniasm -f ${reads} ${paf} 2>/dev/null 1>${gfa}

# Return "good" if the assembly is within a margin.
awk '($1 == "S") { split($4, a, ":"); if (a[3] > 4600000 && a[3] < 5000000) { print "good"} else {print "bad"}}' ${gfa}
