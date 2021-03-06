project(
  'Raptor',
  'cpp',
  version : '0.20.3',
  default_options : [
    'buildtype=release',
    'warning_level=3',
    'cpp_std=c++14'],
  license : 'unlicensed',
  meson_version : '>= 0.48')

############
# CXXFLAGS #
############

raptor_warning_flags = []
cpp = meson.get_compiler('cpp')
foreach cflag: [
  '-Wno-unused-result',
  '-Wno-unused-parameter',
  '-Wno-sign-compare',
  '-Wno-non-virtual-dtor'
  ]
    if cpp.has_argument(cflag)
      raptor_warning_flags += cflag
    endif
endforeach

if cpp.get_id() == 'clang'
  foreach cflag: [
    '-Wno-unused-private-field',
    '-Wno-unused-command-line-argument',
    '-Wno-unused-variable',
    '-Wno-unused-function',
    '-Wno-char-subscripts',
    '-Wno-pessimizing-move',
    '-Wno-format-pedantic',
    '-Wno-keyword-macro'
    ]
      if cpp.has_argument(cflag)
        raptor_warning_flags += cflag
      endif
  endforeach
endif

if cpp.get_id() == 'gcc'
  foreach cflag: [
    '-Wno-unused-variable',
    '-Wno-unused-function',
    '-Wno-char-subscripts',
    '-Wno-pedantic',
    '-Wno-unused-but-set-variable',
    ]
      if cpp.has_argument(cflag)
        raptor_warning_flags += cflag
      endif
  endforeach
endif


raptor_cpp_flags = []
# raptor_cpp_flags = ['-march=native']
# raptor_cpp_flags = ['-march=native', '-DRAPTOR_TESTING_MODE']

opt_testing_mode = get_option('RAPTOR_TESTING_MODE')
if opt_testing_mode
  raptor_cpp_flags += ['-DRAPTOR_TESTING_MODE']
endif

opt_timing_mode = get_option('RAPTOR_DEBUG_TIMINGS')
if opt_timing_mode
  raptor_cpp_flags += ['-DRAPTOR_DEBUG_TIMINGS']
endif

opt_compile_with_pbbam = get_option('WITH_PBBAM')
if opt_compile_with_pbbam
  raptor_cpp_flags += ['-DRAPTOR_COMPILED_WITH_PBBAM']
endif

opt_compile_with_tests = get_option('tests')

################
# Dependencies #
################

# Threads.
raptor_thread_dep = dependency('threads', required : true)

# Zlib.
raptor_zlib_dep = dependency('zlib', required : true)

if opt_compile_with_pbbam
    # Boost. Required for Pbbam.
    boost_dep = dependency('boost', required : true)

    # Pbbam.
    pbbam_dep = dependency('pbbam', fallback : ['pbbam', 'pbbam_dep'])

    # Htslib.
    htslib_dep = dependency('htslib', required : true, version : '>=1.4', fallback : ['htslib', 'htslib_dep'])
else
    boost_dep = []
    pbbam_dep = []
    htslib_dep = []
endif

# Google test.
if (not meson.is_subproject()) and opt_compile_with_tests
  gtest_dep = dependency('gtest', main : true, required : false)

  if not gtest_dep.found()
      gtest_proj = subproject('gtest')
      gtest_inc = gtest_proj.get_variable('gtest_incdir')
      gtest_lib = static_library('gtest', gtest_proj.get_variable('gtest_libsources'),
                                  gtest_proj.get_variable('gtest_mainsources'),
                                  include_directories : gtest_inc)

      gtest_dep = declare_dependency(include_directories : gtest_inc,
                                link_with : gtest_lib, dependencies: raptor_thread_dep)
  endif
endif

#######################
# Configuring headers #
#######################
raptor_version_commit = 'unknown'
git_command = find_program('git', required: false)
if git_command.found()
    git_run = run_command('git', ['log', '-1', '--pretty=%h'])
    if git_run.returncode() == 0
        raptor_version_commit = git_run.stdout().strip()
    endif
endif

raptor_version_h_config = configuration_data()
raptor_version = meson.project_version()
raptor_version_split = meson.project_version().split('.')
raptor_version_h_config.set('RAPTOR_VERSION_MAJOR', raptor_version_split[0])
raptor_version_h_config.set('RAPTOR_VERSION_MINOR', raptor_version_split[1])
raptor_version_h_config.set('RAPTOR_VERSION_PATCH', raptor_version_split[2])
raptor_version_h_config.set('RAPTOR_VERSION_COMMIT', raptor_version_commit)

raptor_version_h = configure_file(
  input : files('include/version.h.in'),
  output : 'version.h',
  configuration : raptor_version_h_config)

###########
# Headers #
###########

raptor_include_directories = [include_directories('.'), include_directories('include'), include_directories('include/ksw2'), include_directories('include/sparsehash')]
# codebase_include_directories = [include_directories('codebase/argparser/src/'), include_directories('codebase/thread_pool/src/'), include_directories('codebase/spoa/include')]
codebase_include_directories = []
codebase_cpp_sources = []
codebase_lib = []

######################
# Sources + codebase #
######################

# subdir('codebase')
subdir('src')

if (not meson.is_subproject()) and opt_compile_with_tests
  subdir('tests')
endif

all_sources = raptor_cpp_sources + codebase_cpp_sources

######################
# The raptor exe.    #
######################

if not meson.is_subproject()
  import('pkgconfig').generate(
    raptor_lib,
    version : meson.project_version(),
    name : 'raptor',
    filebase : 'raptor',
    description : 'Raptor - graph-based mapping of long sequences')
endif


# raptor_dep is not used anywhere?
raptor_dep = declare_dependency(
    include_directories: codebase_include_directories + raptor_include_directories,
    link_with: [raptor_lib, codebase_lib],
    dependencies: [raptor_thread_dep, raptor_zlib_dep, htslib_dep, pbbam_dep, boost_dep],
    version: meson.project_version(),
    compile_args: raptor_warning_flags + raptor_cpp_flags)

if not meson.is_subproject()
  raptor_bin = executable(
    'raptor',
    ['src/raptor/main_raptor.cc'],
    install : true,
    build_by_default : true,
    dependencies : [raptor_thread_dep, raptor_zlib_dep, htslib_dep, pbbam_dep, boost_dep],
    include_directories : codebase_include_directories + raptor_include_directories,
    link_with : [raptor_lib],
    cpp_args : [raptor_warning_flags, raptor_cpp_flags])

  #########################
  # The raptor-index exe. #
  #########################
  # raptor_index_bin = executable(
  #   'raptor-index',
  #   ['src/mains/main_raptor_index.cc'],
  #   dependencies : [raptor_thread_dep, raptor_zlib_dep, boost_dep],
  #   include_directories : codebase_include_directories + raptor_include_directories,
  #   link_with : [raptor_lib],
  #   cpp_args : [raptor_warning_flags, raptor_cpp_flags])

  #########################
  # The raptor-reshape exe. #
  #########################
  raptor_reshape_bin = executable(
    'raptor-reshape',
    ['src/raptor_reshape/main_raptor_reshape.cc'],
    install : true,
    build_by_default : true,
    dependencies : [raptor_thread_dep, raptor_zlib_dep, htslib_dep, pbbam_dep, boost_dep],
    include_directories : codebase_include_directories + raptor_include_directories,
    link_with : [raptor_lib],
    cpp_args : [raptor_warning_flags, raptor_cpp_flags])

  #########################
  # The raptor-fetch exe. #
  #########################
  raptor_fetch_bin = executable(
    'raptor-fetch',
    ['src/raptor_fetch/main_raptor_fetch.cc'],
    install : true,
    build_by_default : true,
    dependencies : [raptor_thread_dep, raptor_zlib_dep, htslib_dep, pbbam_dep, boost_dep],
    include_directories : codebase_include_directories + raptor_include_directories,
    link_with : [raptor_lib],
    cpp_args : [raptor_warning_flags, raptor_cpp_flags])

  ######################
  # Tests              #
  ######################

  if opt_compile_with_tests
    if (gtest_dep.found())
        tests_bin = executable(
            'tests_raptor',
            tests_cpp_sources,
            dependencies : [raptor_thread_dep, raptor_zlib_dep, gtest_dep, htslib_dep, pbbam_dep, boost_dep],
            include_directories : codebase_include_directories + raptor_include_directories + tests_include_directories,
            link_with : [raptor_lib, codebase_lib],
            cpp_args : [raptor_warning_flags, raptor_cpp_flags])
    endif
  endif

endif

# e = executable('testprog', 'test.cc', dependencies : gtest_dep)
# test('gtest test', e)
