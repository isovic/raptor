#! /usr/bin/env python2.7

"""
Example usage:
scripts/graphsim.py test_data/graph/graph-test-6-falcon_assembly_graph/all_ctg.fa test_data/graph/graph-test-6-falcon_assembly_graph/contigs.gfa2 temp/graphsim/graphsim
"""

import re
import sys
import argparse
import fastqparser
import math
import random
import intervaltree
import copy

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

DEBUG_MODE = False
DEBUG_VERBOSE = False

def revcmp(seq):
    return ''.join([COMPLEMENT[v] for v in seq[::-1]])

class GFAGraph:
    def __init__(self, header=None, nodes=None, edges=None):
        self.header = header
        self.nodes = nodes
        self.edges = edges

def load_gfa2(fp_in, add_symmetric_arcs):
    """
    Loads a GFA-2 formatted file and returns the created GFAGraph object.
    Example GFA-2:
        H	VN:Z:2.0
        S	gi|545778205|gb|U00096.3|	4641652	*
        E	edge1	gi|545778205|gb|U00096.3|+	gi|545778205|gb|U00096.3|+	4641652$	4641652$	0	0	*
    """

    def get_coord(string_val):
        if len(string_val) == 0: return 0
        return int(string_val.split('$')[0])
    def parse_header(split_line):
        assert(split_line[0] == 'H')
        header = {}
        for tag in split_line[1:]:
            tag_name, tag_type, tag_val = tag.split(':')
            header[tag_name] = (tag_type, tag_val)
        return header

    def parse_segment(split_line):
        assert(split_line[0] == 'S')
        node = {}
        node['name'] = split_line[1]
        node['len'] = int(split_line[2])
        node['seq'] = split_line[3]
        node['tags'] = {}           # TODO: Parse these.
        node['labels'] = {}
        return node

    def parse_edge(split_line):
        assert(split_line[0] == 'E')
        edge = {}
        edge['name'] = split_line[1]
        edge['v'], edge['v_orient'] = split_line[2][0:-1], split_line[2][-1]
        edge['w'], edge['w_orient'] = split_line[3][0:-1], split_line[3][-1]
        edge['v_start'] = get_coord(split_line[4])
        edge['v_end'] = get_coord(split_line[5])
        edge['w_start'] = get_coord(split_line[6])
        edge['w_end'] = get_coord(split_line[7])
        edge['cigar'] = split_line[8]
        edge['tags'] = {}           # TODO: Parse these.
        edge['labels'] = {}
        return edge

    def create_symmetric_edge(nodes, edge):
        edge_rev = copy.copy(edge)
        edge_rev['name'] += '_rev'

        edge_rev['v'] = edge['w']
        orig_w_len = nodes[edge['w']]['len']
        edge_rev['v_orient'] = '+' if edge['w_orient'] == '-' else '-'
        edge_rev['v_start'] = orig_w_len - edge['w_end']
        edge_rev['v_end'] = orig_w_len - edge['w_start']

        edge_rev['w'] = edge['v']
        orig_v_len = nodes[edge['v']]['len']
        edge_rev['w_orient'] = '+' if edge['v_orient'] == '-' else '-'
        edge_rev['w_start'] = orig_v_len - edge['v_end']
        edge_rev['w_end'] = nodes[edge['v']]['len']- edge['v_start']

        return edge_rev

    header = {}
    nodes = {}
    edges = {}

    for line in fp_in:
        line = line.strip()
        if len(line) == 0: continue
        sl = line.split()
        if sl[0] == 'H':
            current_header_line = parse_header(sl)
            # Check version.
            if 'VN' in current_header_line:
                assert(float(current_header_line['VN'][1]) == 2)
            header.update(current_header_line)
        elif sl[0] == 'S':
            node = parse_segment(sl)
            nodes[node['name']] = node
        elif sl[0] == 'E':
            edge = parse_edge(sl)
            edges[edge['name']] = edge

            if add_symmetric_arcs:
                edge_rev = create_symmetric_edge(nodes, edge)
                edges[edge_rev['name']] = edge_rev

    return GFAGraph(header=header, nodes=nodes, edges=edges)

def make_interval_trees(graph):
    out_edges = {}
    in_edges = {}
    for edge_name, edge in graph.edges.iteritems():
        out_edges.setdefault(edge['v'], [])
        out_edges[edge['v']].append(edge)
        in_edges.setdefault(edge['w'], [])
        in_edges[edge['w']].append(edge)

    interval_dict = {}
    for v, edges in out_edges.iteritems():
        interval_dict[v] = [intervaltree.Interval(edge['v_start'], edge['v_start'] + 1, edge) for edge in edges]

    interval_trees = {}
    for v, intervals in interval_dict.iteritems():
        interval_trees[v] = intervaltree.IntervalTree(intervals)

    return interval_trees

def get_strand_interval(seq, strand, start_pos, end_pos):
    if strand == '+':
        return start_pos, end_pos
    seq_len = len(seq)
    return seq_len - end_pos, seq_len - start_pos

def propagate_path(interval_trees, ref_seqs, ref_seqs_rev, seq_name, seq_strand, pos, remaining_len):
    """
    For a given position on a seq, checks if there is an out edge (or more than one).
    Randomly chooses between all out edges, including the continuation on the current
    sequence.
    After the choice, the position is moved along the chosen sequence, up until either:
        1. Another fork is found.
        2. The end of the sequence is reached and there were no other forks.
        3. The max_dist is reached.
    """

    if remaining_len <= 0:
        if DEBUG_VERBOSE:
            print '=== Stopping (1) ==='
        return None

    # interval_tree = interval_trees[seq_name]
    ref_seq = ref_seqs[seq_name] if seq_strand == '+' else ref_seqs_rev[seq_name]

    def get_fork_pos(interval_trees, seq_name, seq_strand, start, end):
        intervals = interval_trees[seq_name].search(start, end) if seq_name in interval_trees else []
        valid_intervals = [interval.data for interval in intervals if interval.data['v_orient'] == seq_strand]
        valid_intervals = sorted(valid_intervals, key = lambda x: x['v_start'])
        # Find the closest fork position
        fork_pos = end if len(valid_intervals) == 0 else valid_intervals[0]['v_start']
        # Compile a list of outbound edges from this fork position.
        choices = []
        # Find all outbound edges at the closes fork location.
        for i in xrange(0, len(valid_intervals)):
            if valid_intervals[i]['v_start'] != fork_pos:
                break
            choices.append(valid_intervals[i])
        return fork_pos, choices

    fork_pos = None

    if pos >= len(ref_seq):
        # Special case, where there is still some remaining bases to potentially
        # fill, but we reached the end of the sequence.
        # If there are any out edges here, we need to fork out as the only viable
        # option. Otherwise, return None.
        pos = len(ref_seq)  # Intentionally overwriting, to be safe that the coordinate is not bad.
        fork_pos, choices = get_fork_pos(interval_trees, seq_name, seq_strand, pos, pos + 1)
        if DEBUG_VERBOSE:
            print 'Special case: reached the end of the seq.'
    else:
        end_pos = min(len(ref_seq), pos + remaining_len)
        # The +1 below is so that we don't pick already considered forks at the current start
        # position. If there were any forks there, we already made our choice in the previous
        # step, and decided to remain on this sequence.
        fork_pos, choices = get_fork_pos(interval_trees, seq_name, seq_strand, pos + 1, end_pos)
        choices += [None]

    if len(choices) == 0:
        return None

    # Select a direction in which to proceed.
    picked = random.choice(choices)

    if picked == None and pos == len(ref_seq):
        # If we're here, this means we didn't budge from the current sequence,
        # either because we had no out edges, or because the choice was made
        # to stay here.
        if DEBUG_VERBOSE:
            print '=== Stopping (2) ==='
        return None
        # # Either way, the remaining_len should be set to 0, so that the next
        # # iteration will end
        # ret_seq = ref_seq[pos:fork_pos]
        # next_seq_name = seq_name
        # next_seq_strand = seq_strand
        # next_seq_pos = fork_pos
        # next_remaining_len = 0

    elif picked == None:
        ret_seq = ref_seq[pos:fork_pos]
        # Ref coordinates here are on the seq strand instead of fwd to simplify
        # concatenation. This will be converted once all splice segments are generated.
        ret_mapping_paf = ( 'read', len(ret_seq), 0, len(ret_seq),
                            seq_strand, seq_name, len(ref_seq), pos, fork_pos,
                            len(ret_seq), len(ret_seq), 255)
        next_seq_name = seq_name
        next_seq_strand = seq_strand
        next_seq_pos = fork_pos
        next_max_dist = max(0, remaining_len - len(ret_seq))
        if DEBUG_VERBOSE:
            print 'pos = %s, fork_pos = %s, seq_name = %s, seq_len = %d, len(ret_seq) = %d' % (str(pos), str(fork_pos), seq_name, len(ref_seq), len(ret_seq))
            print '=== Staying ==='

    else:
        ret_seq = ref_seq[pos:fork_pos]
        # Ref coordinates here are on the seq strand instead of fwd to simplify
        # concatenation. This will be converted once all splice segments are generated.
        ret_mapping_paf = ( 'read', len(ret_seq), 0, len(ret_seq),
                            seq_strand, seq_name, len(ref_seq), pos, fork_pos,
                            len(ret_seq), len(ret_seq), 255)

        next_seq_name = picked['w']
        next_seq_strand = picked['w_orient']
        next_seq_pos = picked['w_start']
        next_max_dist = max(0, remaining_len - len(ret_seq))
        if DEBUG_VERBOSE:
            print 'pos = %s, fork_pos = %s, seq_name = %s, seq_len = %d, len(ret_seq) = %d' % (str(pos), str(fork_pos), seq_name, len(ref_seq), len(ret_seq))
            print '=== Forking ==='

    return ret_seq, ret_mapping_paf, next_seq_name, next_seq_strand, next_seq_pos, next_max_dist

def generate_exact_insert(trees, ref_seqs, ref_seqs_rev, seq_name, seq_strand, start_pos, sim_read_len, read_prefix='graphsim', zmw_id=0, subread_start=0):
    """
    Extracts a sequence from the reference, as-is, with no mutations added.
    The sequence is composed of spliced parts, extracted from the walk down the graph.
    Graph traversal starts from the sequence and the position specified in this functions
    parameters, and choosing an edge is performed randomly when such an edge is encountered.
    Returns a tuple of: (read_name, read_seq, mappings)
    """

    read_name = 'read-1'
    read_seq = ''

    next_seq_name, next_seq_strand, next_seq_pos, next_remaining_len = seq_name, seq_strand, start_pos + 0, sim_read_len + 0
    mappings = []
    i = 0
    while True:
        if DEBUG_VERBOSE:
            print ''
            print 'Before: next_seq_name = %s, next_seq_strand = %s, next_seq_pos = %d, next_remaining_len = %d' % (next_seq_name, next_seq_strand, next_seq_pos, next_remaining_len)

        new_splice = propagate_path(trees, ref_seqs, ref_seqs_rev, next_seq_name, next_seq_strand, next_seq_pos, next_remaining_len)

        if new_splice == None:
            break

        ret_seq, ret_mapping_paf, next_seq_name, next_seq_strand, next_seq_pos, next_remaining_len = new_splice

        if len(ret_seq) == 0:
            continue
            # Format of ret_mapping = (seq_name, seq_strand, len(ref_seq), pos, fork_pos)
            # Format of ret_mapping_paf:
            #   ret_mapping_paf = ( 'read', len(ret_seq), 0, len(ret_seq),
            #           seq_strand, seq_name, len(ref_seq), pos, fork_pos,
            #           len(ret_seq), len(ret_seq), 255)

        read_seq += ret_seq

        # Check if the new splice simply extends the previous one.
        if len(mappings) > 0 and mappings[-1][4] == ret_mapping_paf[4] and mappings[-1][5] == ret_mapping_paf[5] and mappings[-1][8] == ret_mapping_paf[7]:
            mappings[-1][8] = ret_mapping_paf[8]
        else:
            mappings.append([val for val in ret_mapping_paf])
            mappings[-1][0] = read_name

        # print '  - len(ret_seq) = %d, next_seq_name = %s, next_seq_strand = %s, next_seq_pos = %d, next_remaining_len = %d' % (len(ret_seq), next_seq_name, next_seq_strand, next_seq_pos, next_remaining_len)
        # i += 1
        # if i > 5:
        #     exit(1)
            # break
        # print new_splice

    # Adjust the read mapping coordinates.
    read_name = '{}/{}/{}_{}'.format(read_prefix, zmw_id, subread_start, subread_start + len(read_seq))
    qpos = 0
    for i in xrange(len(mappings)):
        m = mappings[i]
        qspan = m[8] - m[7]
        mappings[i][0] = read_name
        mappings[i][1] = len(read_seq)
        mappings[i][2] = qpos
        mappings[i][3] = qpos + qspan
        # Reverse complement coords should be on the fwd strand in PAF.
        if m[4] == '-':
            mappings[i][7], mappings[i][8] = (m[6] - m[8]), (m[6] - m[7])
        qpos += qspan

    if DEBUG_MODE:
        exit(1)

    return read_name, read_seq, mappings

def generate_mutations(insert_name, insert_seq, insert_mappings,
                        error_rate, frac_snp, frac_ins, frac_del):

    def pick_alt_base(ref_base):
        bases = [val for val in ['A', 'C', 'T', 'G'] if val != ref_base]
        return random.choice(bases)

    def add_op(cigar_list, new_op):
        if len(cigar_list) > 0 and cigar_list[-1][1] == new_op:
            cigar_list[-1][0] += 1
        else:
            cigar_list.append([1, new_op])

    # Define thresholds for random selection of errors.
    thr_snp, thr_ins, thr_del = 0.0, frac_snp, frac_snp + frac_ins

    # Placeholders for return results.
    new_seq = []
    all_cigars = []
    new_mappings = copy.deepcopy(insert_mappings)
    num_eq, num_x, num_ins, num_del = 0, 0, 0, 0

    for m in new_mappings:
        qstart = m[2] + 0
        qend = m[3] + 0

        new_qstart = len(new_seq)
        curr_cigar = []

        i = qstart
        while i < qend:
            ref_base = insert_seq[i]
            # Should we introduce an error here?
            chance = random.random()

            # No error, just continue.
            if chance >= error_rate:
                add_op(curr_cigar, '=')
                new_seq.append(ref_base)
                num_eq += 1
                i += 1
                continue

            # Ok, we need to make an error.
            # Which error do we make?
            chance = random.random()

            if chance >= 0 and chance < thr_ins:
                add_op(curr_cigar, 'X')
                alt_base = pick_alt_base(ref_base)
                new_seq.append(alt_base)
                num_x += 1

            elif chance >= thr_ins and chance < thr_del:
                add_op(curr_cigar, 'I')
                alt_base = pick_alt_base('.')   # Any non ACTG char.
                new_seq.append(alt_base)
                num_ins += 1
                i -= 1

            elif chance >= thr_del and chance < 1.0:
                add_op(curr_cigar, 'D')
                num_del += 1

            i += 1

        new_qend = len(new_seq)

        m[2] = new_qstart
        m[3] = new_qend
        curr_cigar_string = ''.join(['%d%s' % (val[0], val[1]) for val in curr_cigar])
        m.append('cg:Z:%s' % (curr_cigar_string))
        all_cigars.append(curr_cigar)

    # Concatenate the sequence.
    new_seq = ''.join(new_seq)

    # Update the read length.
    for m in new_mappings:
        m[1] = len(new_seq)

    return new_seq, new_mappings, all_cigars, num_eq, num_x, num_ins, num_del

def update_mappings(insert_name, insert_seq, insert_mappings, cigar):
    return insert_name, insert_seq, insert_mappings

def generate_read(trees, ref_seqs, ref_seqs_rev, seq_name, seq_strand, start_pos, sim_read_len,
                    error_rate, frac_snp, frac_ins, frac_del,
                    missing_adapter_rate, missing_adapter_len_lambda,    # If there was a missing adapter, choose the length of the second pass of the SMRT bell.
                    b_rate, b_lambda,
                    read_prefix='graphsim', zmw_id=0, subread_start=0):

                    # missing_adapter_rate, second_pass_mean_len, second_pass_std,    # If there was a missing adapter, choose the length of the second pass of the SMRT bell.

    # This part extracts an exact insert sequence; for example, the fragment which would be part of the SMRT bell.
    insert_name, insert_seq, insert_mappings = generate_exact_insert(trees, ref_seqs, ref_seqs_rev, seq_name, seq_strand, start_pos, sim_read_len, read_prefix='graphsim', zmw_id=zmw_id, subread_start=subread_start)

    chance = random.random()
    if chance < missing_adapter_rate:
        # TODO: Implement this part.
        pass

    read_seq, read_mappings, cigar, num_eq, num_x, num_ins, num_del = generate_mutations(insert_name, insert_seq, insert_mappings,
                                error_rate, frac_snp, frac_ins, frac_del)

    # print 'Cigar:'
    # print ''.join(['%d%s' % (val[0], val[1]) for val in cigar])
    # print 'num_eq = %d, num_x = %d, num_ins = %d, num_del = %d' % (num_eq, num_x, num_ins, num_del)

    # TODO: Left-align indels.

    # read_name, read_seq, read_mappings = update_mappings(insert_name, insert_seq, insert_mappings, cigar)

    print 'Insert_mappings:'
    for m in insert_mappings:
        print m

    return insert_name, read_seq, read_mappings

def mutate_seq():
    mutated = []
    return mutated

def run(ref, gfa, out_prefix, cov, len_mean, len_std, len_min, len_max):
    """
    test_intervals = [intervaltree.Interval(1, 7, 0), intervaltree.Interval(7, 11, 1)]
    test_tree = intervaltree.IntervalTree(test_intervals)
    print 'Test 1 - the end coordinate in search is not inclusive: ', test_tree.search(0, 7)
    print 'Test 2 - the end coordinate in the interval construction is not inclusive: ', test_tree.search(7, 9)
    # Output:
    # Test 1 - the end coordinate in search is not inclusive:  set([Interval(1, 7, 0)])
    # Test 2 - the end coordinate in the interval construction is not inclusive:  set([Interval(7, 11, 1)])
    exit(1)
    """

    # TODO: Parametrize this via the command line:
    error_rate, frac_snp, frac_ins, frac_del = 0.15, 0.25, 0.50, 0.25
    # error_rate, frac_snp, frac_ins, frac_del = 0.00, 0.25, 0.50, 0.25
    missing_adapter_rate, missing_adapter_len_lambda = 0.01, 1.0
    b_rate, b_lambda = 0.0, 1.0

    assert((frac_snp + frac_ins + frac_del) == 1.0)
    assert(len_min > 0)
    assert(len_max > len_min)
    assert(len_max > len_mean)



    ref_seqs = {seq[0][1:].split()[0]: seq[1] for seq in fastqparser.yield_seq([ref])}
    ref_seqs_rev = {key: revcmp(val) for key, val in ref_seqs.iteritems()}

    if not os.path.exists(os.path.dirname(out_prefix)):
        os.makedirs(os.path.dirname(out_prefix))

    if gfa and os.path.exists(gfa):
        with open(gfa, 'r') as fp_in:
            graph = load_gfa2(fp_in, True)
    else:
        graph = GFAGraph(header={}, nodes={}, edges={})

    # print graph.header
    # print graph.nodes
    # print graph.edges

    trees = make_interval_trees(graph)

    total_len = sum([len(seq) for qname, seq in ref_seqs.iteritems()])
    num_generated_reads = 0

    with open('%s.paf' % (out_prefix), 'w') as fp_out_paf, \
            open('%s.fasta' % (out_prefix), 'w') as fp_out_fasta:

        for seq_name, seq in ref_seqs.iteritems():
            total_bases = 0
            target_bases = len(seq) * cov
            seq_len = len(seq)
            while total_bases < target_bases:
                mid_pos = random.randint(0, len(seq))
                seq_strand = random.choice(['+', '-'])

                if DEBUG_MODE:
                    mid_pos = 19188
                    seq_strand = '+'

                sim_read_len = 0
                while sim_read_len < len_min or sim_read_len > len_max:
                    sim_read_len = int(math.floor(random.gauss(len_mean, len_std)))

                total_bases += sim_read_len
                start_pos = max(0, mid_pos - sim_read_len / 2)
                end_pos = min(seq_len, start_pos + sim_read_len)

                if DEBUG_MODE:
                    start = 0
                    mid = 19188
                    end = 100581
                    sim_read_len = 100581
                    total_bases = 100581
                    seq_strand = '+'

                if DEBUG_VERBOSE:
                    print 'start = %d, mid = %d, end = %d, total_bases = %d, sim_read_len = %d, len(seq) = %d, target_bases = %d' % (start_pos, mid_pos, end_pos, total_bases, sim_read_len, len(seq), target_bases)

                new_read = generate_read(trees, ref_seqs, ref_seqs_rev, seq_name, seq_strand, start_pos, sim_read_len,
                                            error_rate, frac_snp, frac_ins, frac_del,
                                            missing_adapter_rate, missing_adapter_len_lambda,
                                            b_rate, b_lambda,
                                            read_prefix='graphsim',
                                            zmw_id=num_generated_reads,
                                            subread_start=0
                                            )

                if new_read == None:
                    continue
                read_name, read_seq, read_mappings = new_read

                fp_out_fasta.write('>%s\n' % (read_name))
                fp_out_fasta.write(read_seq)
                fp_out_fasta.write('\n')

                for m in read_mappings:
                    fp_out_paf.write('\t'.join([str(val) for val in m]))
                    fp_out_paf.write('\n')

                # Just debug output.
                print 'Read mappings:'
                for m in read_mappings:
                    print m[0:12]

                num_generated_reads += 1

                print ''

                # path = generate_path(trees[seq_name], seq_name, seq, start_pos, read_len)
                # propagate_path(trees[seq_name], seq_name, seq, seq_strand, start_pos, read_len)

def parse_args(argv):
    parser = argparse.ArgumentParser(description='Simulates reads from a genome graph, specified by a FASTA/FASTQ file of sequences and a relation graph in GFA-2 format.',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('ref', help='Reference sequences in FASTA/FASTQ format.')
    parser.add_argument('gfa', help='Reference graph in GFA-2 format.')
    parser.add_argument('out_prefix', type=str, default='out.graphsim', help='Prefix of the output files which will be generated.')
    parser.add_argument('--cov', type=float, default=5, help='Coverage of the genome to generate.')
    parser.add_argument('--len-mean', type=float, default=10000, help='Mean read length to simulate.')
    parser.add_argument('--len-std', type=float, default=2000, help='Mean read length to simulate.')
    parser.add_argument('--len-min', type=float, default=500, help='Mean read length to simulate.')
    parser.add_argument('--len-max', type=float, default=70000, help='Mean read length to simulate.')
    # parser.add_argument('--match-rate-mean', type=float, default=, help='Mean read length to simulate.')

    args = parser.parse_args(argv[1:])

    return args

def main(argv=sys.argv):
    args = parse_args(argv)

    run(**vars(args))

if __name__ == "__main__":  # pragma: no cover
    main(sys.argv)          # pragma: no cover
