.PHONY: all clean testing time time2 debug debug-gcc6 data cram unit cram-local cram-external cram-integration tests tools dist configure install build

default: release

clean:
	rm -rf ${BDIR}
	# git clean -xdf .

# Custom define flags:
#	EXPERIMENTAL_QUERY_MASK
#	RAPTOR_TESTING_MODE
#	RAPTOR_DEBUG_TIMINGS

# Default install and build directories:
PREFIX?=${CURDIR}/install

# In case we have a problem with BIN_DIR in cram,
PATH:=${PREFIX}/bin:${PATH}
LD_LIBRARY_PATH:=${PREFIX}/lib64:${PREFIX}/lib:${LD_LIBRARY_PATH}
export PATH LD_LIBRARY_PATH

MESON_FLAGS?="--prefix=${PREFIX} --buildtype=release -DRAPTOR_TESTING_MODE=false -Dc_args=-O3"
BDIR?=meson-release

# Most rules will create BDIR only if it does not already exist.
# ("|" means  "order-only" rule, useful for directory creation.)

build: | ${BDIR}
	ninja -C ${BDIR} reconfigure
	ninja -C ${BDIR}
install: | ${BDIR}
	ninja -C ${BDIR} reconfigure
	ninja -C ${BDIR} install

# This is the only rule that uses MESON_FLAGS.
# If you want to recreate a directory, you can run "make configure", or simply "rm -rf meson-dir".
configure:
	rm -rf ${BDIR} && mkdir -p ${BDIR} && meson ${MESON_FLAGS} ${BDIR}

# These are rules to build specific directories.
# For convenience, you can set "BDIR" to one of these in your shell.
meson-release:
	${MAKE} configure BDIR=$@ \
		MESON_FLAGS="--prefix=${PREFIX} --buildtype=release -DRAPTOR_TESTING_MODE=false -Dc_args=-O3"
meson-testing:
	${MAKE} configure BDIR=$@ \
		MESON_FLAGS="--prefix=${PREFIX} --buildtype=release -DRAPTOR_TESTING_MODE=true -Dc_args=-O3"
meson-time:
	${MAKE} configure BDIR=$@ \
		MESON_FLAGS="--prefix=${PREFIX} --buildtype=release -DRAPTOR_TESTING_MODE=true -DRAPTOR_DEBUG_TIMINGS=true -Dc_args=-O3"
meson-time2:
	${MAKE} configure BDIR=$@ \
		MESON_FLAGS="--prefix=${PREFIX} --buildtype=release -DRAPTOR_DEBUG_TIMINGS=true -Dc_args=-O3"
meson-debug:
	${MAKE} configure BDIR=$@ \
		MESON_FLAGS="--prefix=${PREFIX} --buildtype=debug -Db_sanitize=address -Dc_args=-O3"
meson-debug-gcc6:
	${MAKE} configure BDIR=$@ \
		MESON_FLAGS="--prefix=${PREFIX} --buildtype=debug -Db_sanitize=address -Dc_args=-O3"

# These rules ignore your current BDIR setting, but they rely on $PREFIX via $MESON_FLAGS.
# They all reconfigure, rebuild, and install.
release: | meson-release
	${MAKE} install BDIR=meson-release
testing: | meson-testing
	${MAKE} install BDIR=meson-testing
time: | meson-time
	${MAKE} install BDIR=meson-time
time2: | meson-time2
	${MAKE} install BDIR=meson-time
debug: | meson-debug
	${MAKE} install BDIR=meson-debug
debug-gcc6: | meson-debug-gcc6
	${MAKE} install BDIR=meson-debug-gcc6

dist: release
	cd meson-release && ninja-dist

###########################################
### Tests.                              ###
###########################################
third-party/cram/scripts/cram:
	git submodule update --init third-party/cram

data: raptor-test-data/README.md
	cd raptor-test-data && git pull

raptor-test-data/README.md:
	git clone https://github.com/isovic/raptor-test-data.git

cram: installed third-party/cram/scripts/cram cram-local #cram-external

unit: release
	meson-release/tests_raptor

unit-testing: testing
	meson-testing/tests_raptor

cram-local: installed
	scripts/cram -E tests/cram/local/*.t tests/cram/local-graph/*.t

cram-external: installed raptor-test-data/README.md
	scripts/cram -E tests/cram/external/*.t

cram-integration: installed tools/miniasm/miniasm raptor-test-data/README.md
	scripts/cram -E tests/cram/integration/*.t

installed: install/bin/raptor  # see BIN_DIR in scripts/cram

tests: unit cram



###########################################
### Tools. Required only for testing.   ###
###########################################
tools: tools/miniasm/miniasm

tools/miniasm/miniasm:
	mkdir -p tools && cd tools && git clone https://github.com/lh3/miniasm && cd miniasm && make
