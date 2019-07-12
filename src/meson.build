###########
# sources #
###########

raptor_cpp_sources = files([
  'aligner/aligner_edlib.cc',
  'aligner/aligner_factory.cc',
  'aligner/aligner_ksw2_double.cc',
  'aligner/aligner_ksw2_single.cc',
  'aligner/aligner_util.cc',
  'aligner/alignment_result.cc',
  'aligner/anchor_aligner.cc',
  'aligner/cigar.cc',
  'aligner/sam_parser.cc',
  'aligner/difflib_edlib.cc',

#  'containers/chain.cc',
#  'containers/chained_alignment.cc',
  'containers/region/region_aligned.cc',
  'containers/region/region_mapped.cc',
  'containers/mapping_result/linear_mapping_result.cc',
  'containers/mapping_result/graph_mapping_result.cc',
  'containers/mapping_result/aligned_mapping_result.cc',
  'containers/mapping_env.cc',
  'containers/path_alignment.cc',

#  'containers/raptor_chain_alignment_result.cc',
  'debug_tools/write_seed_hit_1.cc',

  'graph/anchor_graph_edge.cc',
  'graph/local_path.cc',
  'graph/segment_edge.cc',
  'graph/segment_graph_parser.cc',
  'graph/segment_graph.cc',
  'graph/segment_node.cc',
  'graph/split_segment_graph.cc',

  'graph_aligner/anchor_graph_aligner.cc',

  'ksw2/kalloc.cc',
  'ksw2/ksw2_extd2_sse.cc',
  'ksw2/ksw2_exts2_sse.cc',
  'ksw2/ksw2_extz2_sse.cc',
  'ksw2/ksw2_ll_sse.cc',

  'lib/edlib.cc',
  'lib/thread_pool.cc',
  'lib/argparser.cc',

  'log/log_tools.cc',
  'log/log_system.cc',

  'index/index_params.cc',
  'index/minimizer_generator.cc',
  'index/minimizer_index_util.cc',
  'index/minimizer_index.cc',

  'params/params_aligner.cc',
  'params/params_mapper.cc',
  'params/params_raptor.cc',
  'params/params_raptor_index.cc',

  'raptor/backtrack_list_node.cc',
#  'raptor/chain_aligner.cc',
  'raptor/dp_chain.cc',
  'raptor/graph_mapper.cc',
  'raptor/graph_mapper_tools.cc',
  'raptor/index_factory.cc',
  'raptor/interval_tree_builder.cc',
  'raptor/mapper.cc',
  'raptor/mapper_tools.cc',
  'raptor/path_aligner.cc',
  'raptor/raptor_aligner.cc',
  'raptor/raptor.cc',

  'raptor_reshape/params_raptor_reshape.cc',

  'raptor_fetch/params_raptor_fetch.cc',
  'raptor_fetch/args_raptor_fetch.cc',
  'raptor_fetch/overlaps/overlap.cc',
  'raptor_fetch/overlaps/overlap_compact.cc',
  'raptor_fetch/overlaps/overlap_file.cc',

  'sequences/sequence_file_parser_factory.cc',
  'sequences/sequence_file_parser_fastx.cc',
  'sequences/sequence_file_parser_gfa1.cc',
  'sequences/sequence_file_parser_gfa2.cc',
  'sequences/sequence_file_parser_sam.cc',
  'sequences/sequence_file_parser_bam.cc',
  'sequences/sequence_file_parser_utils.cc',
  'sequences/random_access_sequence_file.cc',
  'sequences/sequence_file.cc',
  'sequences/sequence.cc',
  'sequences/sequence_serializer.cc',
  'sequences/sequence_file_composite_fofn.cc',
  'sequences/sequence_file_composite_pbxml.cc',

  'utility/memtime.cc',
  'utility/fofn.cc',
  'utility/stringutil.cc',
  'utility/tictoc.cc',

  'writer/output_formatter.cc',
  'writer/raptor_results_writer_factory.cc',
  'writer/raptor_results_writer_stream.cc',
  'writer/raptor_results_writer_bam.cc',

  'raptor/args_raptor.cc',
  'raptor_reshape/args_raptor_reshape.cc',
  'raptor_fetch/args_raptor_fetch.cc',

])

raptor_extra_flags = []

raptor_lib_install = (not meson.is_subproject()) or (get_option('default_library') == 'shared')

raptor_lib = library(
  'raptor',
  raptor_cpp_sources,
  soversion : 0,
  version : meson.project_version(),
  install : raptor_lib_install,
  link_with : codebase_lib,
  dependencies : [raptor_thread_dep, raptor_zlib_dep, htslib_dep, pbbam_dep, boost_dep],
  include_directories : raptor_include_directories + codebase_include_directories,
  cpp_args : [raptor_extra_flags, raptor_warning_flags, raptor_cpp_flags])
