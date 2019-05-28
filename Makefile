.PHONY: all clean testing time time2 debug debug-gcc6 data cram unit cram-local cram-external cram-integration tests tools dist configure install rebuild build

all: build  # for consumers expecting to see a build/ directory, but should be build-release

clean:
	rm -rf ${BDIR}
	# git clean -xdf .

# Custom define flags:
#	EXPERIMENTAL_QUERY_MASK
#	RAPTOR_TESTING_MODE
#	RAPTOR_DEBUG_TIMINGS

PREFIX?=${CURDIR}/PREFIX
BDIR?=build
# Most rules will create BDIR only if it does not already exist.
# ("|" means  "order-only" rule, useful for directory creation.)

rebuild: | ${BDIR}
	ninja -C ${BDIR} reconfigure
	ninja -C ${BDIR}
install: | ${BDIR}
	ninja -C ${BDIR} reconfigure
	ninja -C ${BDIR} install

# "meson --reconfigure" is not idempotent, but
# "ninja reconfigure" works even the first time.

MESON_FLAGS?=--buildtype=debug -DRAPTOR_TESTING_MODE=false -Dc_args=-O0 --prefix=${PREFIX}

# This is the only rule that uses MESON_FLAGS.
# If you want to recreate a directory, you can run "make configure", or simply "rm -rf build-dir".
configure:
	rm -rf ${BDIR} && mkdir -p ${BDIR} && meson ${MESON_FLAGS} ${BDIR}

# These are rules to build specific directories.
# For convenience, you can set "BDIR" to one of these in your shell.
build: # default expected by old ipa/
	${MAKE} configure BDIR=$@
build-release:
	${MAKE} configure BDIR=$@ \
		MESON_FLAGS="--buildtype=release -DRAPTOR_TESTING_MODE=false -Dc_args=-O3"
build-testing:
	${MAKE} configure BDIR=$@ \
		MESON_FLAGS="--buildtype=release -DRAPTOR_TESTING_MODE=true -Dc_args=-O3"
build-time:
	${MAKE} configure BDIR=$@ \
		MESON_FLAGS="--buildtype=release -DRAPTOR_TESTING_MODE=true -DRAPTOR_DEBUG_TIMINGS=true -Dc_args=-O3"
build-time2:
	${MAKE} configure BDIR=$@ \
		MESON_FLAGS="--buildtype=release -DRAPTOR_DEBUG_TIMINGS=true -Dc_args=-O3"
build-debug:
	${MAKE} configure BDIR=$@ \
		MESON_FLAGS="--buildtype=debug -Db_sanitize=address -Dc_args=-O3"
build-debug-gcc6:
	${MAKE} configure BDIR=$@ \
		MESON_FLAGS="--buildtype=debug -Db_sanitize=address -Dc_args=-O3"

# These rules ignore your current BDIR setting.
release: | build-release
	${MAKE} build BDIR=build-release
testing: | build-testing
	${MAKE} build BDIR=build-testing
time: | build-time
	${MAKE} build BDIR=build-time
time2: | build-time2
	${MAKE} build BDIR=build-time
debug: | build-debug
	${MAKE} build BDIR=build-debug
debug-gcc6: | build-debug-gcc6
	${MAKE} build BDIR=build-debug-gcc6

build/raptor: release

build-testing/raptor: testing

dist: release
	cd build && ninja-dist

###########################################
### Tests.                              ###
###########################################
data: raptor-test-data/README.md
	cd raptor-test-data && git pull

raptor-test-data/README.md:
	git clone https://github.com/isovic/raptor-test-data.git

cram: build/raptor third-party/cram/scripts/cram cram-local cram-external

unit: build/raptor
	build/tests_raptor

unit-testing: build-testing/raptor
	build-testing/tests_raptor

cram-local: build/raptor
	scripts/cram tests/cram/local/*.t tests/cram/local-graph/*.t

cram-external: build/raptor raptor-test-data/README.md
	scripts/cram tests/cram/external/*.t

cram-integration: build/raptor tools/miniasm/miniasm raptor-test-data/README.md
	scripts/cram tests/cram/integration/*.t

tests: unit cram



###########################################
### Tools. Required only for testing.   ###
###########################################
tools: tools/miniasm/miniasm

tools/miniasm/miniasm:
	mkdir -p tools && cd tools && git clone https://github.com/lh3/miniasm && cd miniasm && make
