.PHONY: all clean testing time time2 debug debug-gcc6 data cram unit cram-local cram-external cram-integration tests tools dist configure install rebuild

default: | build-release
	${MAKE} rebuild

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

MESON_FLAGS?= --prefix=${PREFIX} --buildtype=debug -DRAPTOR_TESTING_MODE=false -Dc_args=-O0

# This is the only rule that uses MESON_FLAGS.
# If you want to recreate a directory, you can run "make configure", or simply "rm -rf build-dir".
configure:
	rm -rf ${BDIR} && mkdir -p ${BDIR} && meson ${MESON_FLAGS} ${BDIR}

# These are rules to build specific directories.
# For convenience, you can set "BDIR" to one of these in your shell.
build:
	${MAKE} configure BDIR=$@
build-release:
	${MAKE} configure BDIR=$@ \
		MESON_FLAGS="--prefix=${PREFIX} --buildtype=release -DRAPTOR_TESTING_MODE=false -Dc_args=-O3"
build-testing:
	${MAKE} configure BDIR=$@ \
		MESON_FLAGS="--prefix=${PREFIX} --buildtype=release -DRAPTOR_TESTING_MODE=true -Dc_args=-O3"
build-time:
	${MAKE} configure BDIR=$@ \
		MESON_FLAGS="--prefix=${PREFIX} --buildtype=release -DRAPTOR_TESTING_MODE=true -DRAPTOR_DEBUG_TIMINGS=true -Dc_args=-O3"
build-time2:
	${MAKE} configure BDIR=$@ \
		MESON_FLAGS="--prefix=${PREFIX} --buildtype=release -DRAPTOR_DEBUG_TIMINGS=true -Dc_args=-O3"
build-debug:
	${MAKE} configure BDIR=$@ \
		MESON_FLAGS="--prefix=${PREFIX} --buildtype=debug -Db_sanitize=address -Dc_args=-O3"
build-debug-gcc6:
	${MAKE} configure BDIR=$@ \
		MESON_FLAGS="--prefix=${PREFIX} --buildtype=debug -Db_sanitize=address -Dc_args=-O3"

# These rules ignore your current BDIR setting, but they rely on $PREFIX via $MESON_FLAGS.
# They all reconfigure, rebuild, and install.
release: | build-release
	${MAKE} rebuild BDIR=build-release
testing: | build-testing
	${MAKE} rebuild BDIR=build-testing
time: | build-time
	${MAKE} rebuild BDIR=build-time
time2: | build-time2
	${MAKE} rebuild BDIR=build-time
debug: | build-debug
	${MAKE} rebuild BDIR=build-debug
debug-gcc6: | build-debug-gcc6
	${MAKE} rebuild BDIR=build-debug-gcc6

build-testing/raptor: testing

dist: release
	cd build-release && ninja-dist

###########################################
### Tests.                              ###
###########################################
third-party/cram/scripts/cram:
	git submodule update --init third-party/cram

data: raptor-test-data/README.md
	cd raptor-test-data && git pull

raptor-test-data/README.md:
	git clone https://github.com/isovic/raptor-test-data.git

cram: build/raptor third-party/cram/scripts/cram cram-local #cram-external

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
