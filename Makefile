.PHONY: all clean testing time time2 debug debug-gcc6 data cram unit cram-local cram-external cram-integration tests tools dist

all: release

clean:
	rm -rf build

# Custom define flags:
#	EXPERIMENTAL_QUERY_MASK
#	RAPTOR_TESTING_MODE
#	RAPTOR_DEBUG_TIMINGS

release:
	@echo "[Invoking Meson]"
	@(cd build && meson --reconfigure && ninja) || (mkdir -p build && cd build && meson --buildtype=release -DRAPTOR_TESTING_MODE=false -Dc_args=-O3 && ninja)

testing:
	@echo "[Invoking Meson]"
	@(cd build-testing && meson --reconfigure && ninja) || (mkdir -p build-testing && cd build-testing && meson --buildtype=release -DRAPTOR_TESTING_MODE=true -Dc_args=-O3 && ninja)

time:
	@echo "[Invoking Meson]"
	@(cd build-time && meson --reconfigure && ninja) || (mkdir -p build-time && cd build-time && meson --buildtype=release -DRAPTOR_TESTING_MODE=true -DRAPTOR_DEBUG_TIMINGS=true -Dc_args=-O3 && ninja)

time2:
	@echo "[Invoking Meson]"
	@(cd build-time2 && meson --reconfigure && ninja) || (mkdir -p build-time2 && cd build-time2 && meson --buildtype=release -DRAPTOR_DEBUG_TIMINGS=true -Dc_args=-O3 && ninja)

debug:
	@echo "[Invoking Meson]"
	@(cd build-debug && meson --reconfigure && ninja) || (mkdir -p build-debug && cd build-debug && (meson --buildtype=debug -Db_sanitize=address -Dc_args=-O3) && ninja)

debug-gcc6:
	@echo "[Invoking Meson]"
	@(cd build-debug-gcc6 && meson --reconfigure && ninja) || (mkdir -p build-debug-gcc6 && cd build-debug-gcc6 && (env CC=gcc-6 CXX=g++-6 meson --buildtype=debug -Db_sanitize=address -Dc_args=-O3) && ninja)

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
