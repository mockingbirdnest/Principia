ARGS	:= -std=c++23 -stdlib=libc++ -O3 -g --system-header-prefix=serialization/ --system-header-prefix=absl/ --system-header-prefix=benchmark/ --system-header-prefix=google/ --system-header-prefix=gtest/ --system-header-prefix=gmock/ --system-header-prefix=zfp/ -fPIC -fbracket-depth=257 -fdeclspec -fexceptions -ferror-limit=1000 -ffp-contract=off -fno-omit-frame-pointer -fno-char8_t -Wall -Wpedantic -Wno-c++23-extensions -Wno-deprecated-declarations -Wno-deprecated-literal-operator -Wno-char-subscripts -Wno-elaborated-enum-class -Wno-gnu-anonymous-struct -Wno-c99-extensions -Wno-gnu-zero-variadic-macro-arguments -Wno-mathematical-notation-identifier-extension -Wno-nested-anon-types -Wno-nullability-extension -Wno-unknown-pragmas -Wno-unused-command-line-argument -DPROJECT_DIR='std::filesystem::path("ksp_plugin_adapter/")' -DSOLUTION_DIR='std::filesystem::path("./")' -DTEMP_DIR='std::filesystem::path("/tmp")' -DPLATFORM_WITH_CPU_FEATURES='"x64"' -DGTEST_HAS_ABSL -DNDEBUG -DPRINCIPIA_REQUIRES_AVX=0 -DPRINCIPIA_REQUIRES_FMA=0 -m64 -msse4.1 -I. -Ideps/protobuf/src -Ideps/protobuf/third_party/utf8_range -Ideps/gipfeli/include -Ideps/abseil-cpp/install/include -Ideps/re2 -Ideps/core-math/include -Ideps/zfp/include -Ideps/config/include -Ideps/multiprecision/include -Ideps/googletest/googlemock/include -Ideps/googletest/googletest/include -Ideps/googletest/googlemock/ -Ideps/googletest/googletest/ -Ideps/lodepng/.. -Ideps/benchmark/include


PRINCIPIA_MODULE_INTERFACE_UNITS := $(wildcard */*.ixx)
PRINCIPIA_HEADERS := $(wildcard */*.hpp)
PRINCIPIA_TRANSLATION_UNITS := $(wildcard */*.cpp)
PRINCIPIA_TEST_TRANSLATION_UNITS := $(wildcard */*_test.cpp)
PRINCIPIA_TEST_BINS := $(PRINCIPIA_TEST_TRANSLATION_UNITS:%.cpp=bin/%)
PROTO_FILES := $(wildcard */*.proto)
PROTO_TRANSLATION_UNITS                 := $(PROTO_FILES:.proto=.pb.cc)
PROTO_HEADERS                           := $(PROTO_FILES:.proto=.pb.h)
PROTO_PRECOMPILED_HEADERS               := $(PROTO_HEADERS:%=phu/%.pcm)
PROTO_OBJECTS := $(addprefix obj/, $(PROTO_TRANSLATION_UNITS:.cc=.cc.o))
ABSL_HEADERS := $(wildcard deps/abseil-cpp/install/include/absl/*/*.h)
ABSL_PRECOMPILED_HEADERS := $(ABSL_HEADERS:deps/abseil-cpp/install/include/%.h=phu/%.h.pcm)
GIPFELI_HEADERS := $(wildcard deps/gipfeli/include/*/*.h)
GIPFELI_PRECOMPILED_HEADERS := $(GIPFELI_HEADERS:deps/gipfeli/include/%.h=phu/%.h.pcm)
PROTOBUF_HEADERS := $(wildcard deps/protobuf/src/google/protobuf/*.h) $(wildcard deps/protobuf/src/google/protobuf/*/*.h)
PROTOBUF_PRECOMPILED_HEADERS := $(PROTOBUF_HEADERS:deps/protobuf/src/%.h=phu/%.h.pcm)
GMOCK_PRECOMPILED_HEADERS := phu/gmock/gmock.h.pcm
GTEST_PRECOMPILED_HEADERS := phu/gtest/gtest.h.pcm

dependencies: compute_import_dependencies.py \
	$(PRINCIPIA_MODULE_INTERFACE_UNITS) \
	$(PRINCIPIA_HEADERS) \
	$(PRINCIPIA_TRANSLATION_UNITS)
	python3 $< > dependencies

.PHONY: clean

clean:
	rm -rf phu/ # Precompiled header units.
	rm -rf pmu/ # Precompiled module interface units.
	rm -rf obj/ # Objects.
	rm -rf bin/ # Binaries.

# Dependencies

ifneq ($(MAKECMDGOALS), clean)
include dependencies
endif

$(PROTO_OBJECTS): obj/%.o: %

pmu/std.pcm: /usr/lib/llvm-23/share/libc++/v1/std.cppm
	@mkdir -p $(@D)
	clang++-23 $(ARGS) --precompile $< -o $@
	
phu/%.pcm:
	@mkdir -p $(@D)
	clang++-23 $(ARGS) -fmodule-header $< -o $@


# Note the -xc++-module because Clang does not know about the .ixx extension.
pmu/%.pcm:
	@mkdir -p $(@D)
	clang++-23 $(ARGS) \
		$(addprefix -fmodule-file=, $(filter phu/%.pcm, $^)) \
		$(foreach module, $(filter pmu/%.pcm, $^), -fmodule-file=$(module:pmu/%.pcm=%)=$(module)) \
		--precompile -xc++-module $< -o $@

$(PROTO_PRECOMPILED_HEADERS): phu/%.pcm: %
$(ABSL_PRECOMPILED_HEADERS): phu/%.pcm: deps/abseil-cpp/install/include/%
$(GIPFELI_PRECOMPILED_HEADERS): phu/%.pcm: deps/gipfeli/include/%
$(GMOCK_PRECOMPILED_HEADERS): phu/%.pcm: deps/googletest/googlemock/include/%
$(GTEST_PRECOMPILED_HEADERS): phu/%.pcm: deps/googletest/googletest/include/%
$(PROTOBUF_PRECOMPILED_HEADERS): phu/%.pcm: deps/protobuf/src/%

# We don't do dependency resolution on the protos; we compile them all at once.
$(PROTO_HEADERS) $(PROTO_TRANSLATION_UNITS): $(PROTO_FILES)
	deps/protobuf/protoc \
	-I deps/protobuf/src/ -I . $^ --cpp_out=.

obj/gmock_main.o: deps/googletest/googlemock/src/gmock_main.cc
	clang++-23 $(ARGS) --compile $< -o $@

obj/gmock-all.o: deps/googletest/googlemock/src/gmock-all.cc
	clang++-23 $(ARGS) --compile $< -o $@
obj/gtest-all.o: deps/googletest/googletest/src/gtest-all.cc
	clang++-23 $(ARGS) --compile $< -o $@

$(PRINCIPIA_TEST_BINS): obj/gmock_main.o obj/gmock-all.o obj/gtest-all.o $(PROTO_OBJECTS)

obj/%.o: 
	@mkdir -p $(@D)
	clang++-23 $(ARGS) \
		$(addprefix -fmodule-file=, $(filter phu/%.pcm, $^)) \
		$(foreach module, $(filter pmu/%.pcm, $^), -fmodule-file=$(module:pmu/%.pcm=%)=$(module)) \
		--compile $< -o $@

bin/%:
	@mkdir -p $(@D)
	clang++-23 $(ARGS) $^ \
		deps/protobuf/libprotobuf.a \
		deps/protobuf/third_party/utf8_range/libutf8*.a \
		deps/gipfeli/libgipfeli.a \
		-Wl,--start-group deps/abseil-cpp/install/lib/libabsl_*.a -Wl,--end-group \
		deps/re2/libre2.a \
		-o $@
