#! /usr/bin/env python3

import sys
import argparse
import collections
import intervaltree

MappingTuple = collections.namedtuple('MappingTuple', ['qname', 'qrev', 'qstart', 'qend', 'qlen', 'tname', 'trev', 'tstart', 'tend', 'tlen'])

RESULT_TYPE_NONE = 'none'
RESULT_TYPE_MISS = 'miss'
RESULT_TYPE_FULL = 'full'
RESULT_TYPE_PARTIAL = 'partial'

################################
######### Utility tools ########
import sys
import os
from time import gmtime, strftime
import traceback
import fnmatch
### Logs messages to STDERR and an output log file if provided (opened elsewhere).
def log(message, fp_log=sys.stderr):
  timestamp = strftime("%Y/%m/%d %H:%M:%S", gmtime());      # pragma: no cover
  if (fp_log != None):                                      # pragma: no cover
      fp_log.write('[%s] %s\n' % (timestamp, message))      # pragma: no cover
      fp_log.flush();                                       # pragma: no cover
################################

COMPLEMENT = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', 'a': 't', 'c': 'g', 't': 'a', 'g': 'c'}

def hash_mappings(truth_fn):
    mappings_dict = collections.defaultdict(list)
    with open(truth_fn) as fp_in:
        for line in fp_in:
            sl = line.strip().split()
            mapping = MappingTuple(sl[0], False, int(sl[2]), int(sl[3]), int(sl[1]),
                                    sl[5], True if sl[4] == '-' else False, int(sl[7]), int(sl[8]), int(sl[6]))
            mappings_dict[mapping.qname].append(mapping)
    return mappings_dict

def calc_overlap_span(s1, e1, s2, e2):
    if s2 >= e1 or e2 <= s1:
        return 0
    os = max(s1, s2)
    oe = min(e1, e2)
    return oe - os

def evaluate_mappings(qname, q_truth_mappings, q_test_mappings, min_ovl_frac):
    """
    Algorithm:
        - Input all truth mappings for a single qname, and all test mappings for the same qname.
        - Build interval trees for the truth mappings.
        - Define an empty fraction accumulator for each truth mapping (i.e. [0] * len(q_truth_mappings) ).
        - For each test mapping:
            - Find overlaps with the truth mappings in the interval tree.
            - Compute the overlapping span between the truth and the test mapping in both query and target coords.
            - Compute fraction of the mapping span covered by the test mapping.
            - Take the minimum fraction of the query or target coordinates.
            - Add the selected fraction to an accumulator of this particular truth mapping.
        - Count the accumulator bins >= min_ovl_frac treshold.
        - If all bins have value >= min_ovl_frac, then report a fully mapped query.
        - If 0 bins have value >= min_ovl_frac, then report a missed mapping.
        - Otherwise, it's a partial mapping.

    """

    interval_dict = collections.defaultdict(list)
    for i, mapping in enumerate(q_truth_mappings):
        target = mapping.tname + ('-' if mapping.trev else '+')
        interval_dict[target].append(intervaltree.Interval(mapping.tstart, mapping.tend, i))

    interval_trees = {}
    for target, intervals in interval_dict.items():
        interval_trees[target] = intervaltree.IntervalTree(intervals)

    result_array = [0] * len(q_truth_mappings)
    result_array_hit_count = [0] * len(q_truth_mappings)

    # Inspect each query portion separately.
    for test_mapping in q_test_mappings:

        # Find all hits for this query.
        target = test_mapping.tname + ('-' if test_mapping.trev else '+')
        if target not in interval_trees:
            continue
        overlaps = interval_trees[target].overlap(test_mapping.tstart, test_mapping.tend)

        # Check the fraction of each hit, and accumulate.
        for hit in overlaps:
            truth_id = hit.data
            truth_mapping = q_truth_mappings[truth_id]

            # Calculate the amount of overlap between the query and target coordinates.
            target_overlap_span = calc_overlap_span(truth_mapping.tstart, truth_mapping.tend, test_mapping.tstart, test_mapping.tend)
            query_overlap_span = calc_overlap_span(truth_mapping.qstart, truth_mapping.qend, test_mapping.qstart, test_mapping.qend)

            target_overlap_frac = min(float(target_overlap_span) / float(truth_mapping.tend - truth_mapping.tstart), float(target_overlap_span) / float(test_mapping.tend - test_mapping.tstart))
            query_overlap_frac = min(float(query_overlap_span) / float(truth_mapping.qend - truth_mapping.qstart), float(query_overlap_span) / float(test_mapping.qend - test_mapping.qstart))

            result_array[truth_id] += min(target_overlap_frac, query_overlap_frac)
            result_array_hit_count[truth_id] += 1

    num_hits = 0
    for val in result_array:
        if val >= min_ovl_frac:
            num_hits += 1

    type_ = RESULT_TYPE_NONE
    if num_hits == 0 and (len(q_truth_mappings) == 0 or len(q_test_mappings) == 0):
        type_ = RESULT_TYPE_NONE
    elif num_hits == 0:
        type_ = RESULT_TYPE_MISS
    elif num_hits == len(result_array):
        type_ = RESULT_TYPE_FULL
    else:
        type_ = RESULT_TYPE_PARTIAL

    return type_, result_array, result_array_hit_count

def run(truth_fn, test_fn_list, min_ovl_frac):
    all_truth_mappings = hash_mappings(truth_fn)
    all_results = []

    for test_fn in test_fn_list:
        all_test_mappings = hash_mappings(test_fn)

        set_missing_in_truth = set()
        set_missing_in_test = set()

        result_count = collections.defaultdict(int)

        for qname, q_test_mappings in all_test_mappings.items():
            expected = all_truth_mappings[qname]
            if not expected:
                log('Warning: no mappings in the truth set for qname: "{}". Skipping.\n'.format(qname))
                set_missing_in_truth.add(qname)

        for qname, q_truth_mappings in all_truth_mappings.items():
            q_test_mappings = all_test_mappings[qname]
            if not q_test_mappings:
                set_missing_in_test.add(qname)
                continue

            result_type, result_array, result_array_hit_count = evaluate_mappings(qname, q_truth_mappings, q_test_mappings, min_ovl_frac)
            result_count[result_type] += 1

        curr_results = [test_fn, len(all_truth_mappings)] + [result_count[type_] for type_ in [RESULT_TYPE_FULL, RESULT_TYPE_PARTIAL, RESULT_TYPE_MISS, RESULT_TYPE_NONE]] + [len(set_missing_in_truth), len(set_missing_in_test)]
        all_results.append(curr_results)

    header = ['test_fn', 'truth', 'full', 'partial', 'miss', 'none', 'missing_in_truth', 'missing_in_test']
    sys.stdout.write('{}\n'.format('\t'.join(header)))
    for vals in all_results:
        sys.stdout.write('{}\n'.format('\t'.join([str(val) for val in vals])))

def parse_args(argv):
    parser = argparse.ArgumentParser(description='Compares a test set of mappings to a truth set.',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--min-ovl-frac', type=float, default=0.10, help='Minimum fraction overlap between the test overlap and the truth overlap.')
    parser.add_argument('truth_fn', help='Input file (paf) containing "true" mappings.')
    parser.add_argument('test_fn_list', nargs='*', help='A set of mappings to be tested')

    args = parser.parse_args(argv[1:])

    return args

def main(argv=sys.argv):
    args = parse_args(argv)

    run(**vars(args))

if __name__ == "__main__":  # pragma: no cover
    main(sys.argv)          # pragma: no cover
