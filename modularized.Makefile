ARGS	:= -std=c++23 -stdlib=libc++ -O3 -g --system-header-prefix=serialization/ --system-header-prefix=absl/ --system-header-prefix=benchmark/ --system-header-prefix=google/ --system-header-prefix=gtest/ --system-header-prefix=gmock/ --system-header-prefix=zfp/ -fPIC -fbracket-depth=257 -fdeclspec -fexceptions -ferror-limit=1000 -ffp-contract=off -fno-omit-frame-pointer -fno-char8_t -Wall -Wpedantic -Wno-c++23-extensions -Wno-deprecated-declarations -Wno-deprecated-literal-operator -Wno-char-subscripts -Wno-elaborated-enum-class -Wno-gnu-anonymous-struct -Wno-c99-extensions -Wno-gnu-zero-variadic-macro-arguments -Wno-mathematical-notation-identifier-extension -Wno-nested-anon-types -Wno-nullability-extension -Wno-unknown-pragmas -Wno-unused-command-line-argument -DPROJECT_DIR='std::filesystem::path("ksp_plugin_adapter/")' -DSOLUTION_DIR='std::filesystem::path("./")' -DTEMP_DIR='std::filesystem::path("/tmp")' -DPLATFORM_WITH_CPU_FEATURES='"x64"' -DGTEST_HAS_ABSL -DNDEBUG -DPRINCIPIA_REQUIRES_AVX=0 -DPRINCIPIA_REQUIRES_FMA=0 -m64 -msse4.1 -I. -Ideps/protobuf/src -Ideps/protobuf/third_party/utf8_range -Ideps/gipfeli/include -Ideps/abseil-cpp/install/include -Ideps/re2 -Ideps/core-math/include -Ideps/zfp/include -Ideps/config/include -Ideps/multiprecision/include -Ideps/googletest/googlemock/include -Ideps/googletest/googletest/include -Ideps/googletest/googlemock/ -Ideps/googletest/googletest/ -Ideps/lodepng/.. -Ideps/benchmark/include

clean:
	rm -rf phu/ # Precompiled header units.
	rm -rf pmu/ # Precompiled module units.
	rm -rf obj/ # Objects.
	rm -rf bin/ # Binaries.

PRINCIPIA_MODULE_INTERFACE_UNITS := $(wildcard */*.ixx)
PRINCIPIA_PRECOMPILED_HEADERS := phu/base/macros.pcm 
SERIALIZATION_HEADERS := $(wildcard serialization/*.pb.h)
SERIALIZATION_PRECOMPILED_HEADERS := $(addprefix phu/, $(SERIALIZATION_HEADERS:.h=.pcm))
SERIALIZATION_TRANSLATION_UNITS := $(wildcard serialization/*.pb.cc)
SERIALIZATION_OBJECTS := $(addprefix obj/, $(SERIALIZATION_TRANSLATION_UNITS:.cc=.o))
ABSL_HEADERS := $(wildcard deps/abseil-cpp/install/include/absl/*/*.h)
ABSL_PRECOMPILED_HEADERS := $(ABSL_HEADERS:deps/abseil-cpp/install/include/%.h=phu/%.pcm)
GIPFELI_HEADERS := $(wildcard deps/gipfeli/include/*/*.h)
GIPFELI_PRECOMPILED_HEADERS := $(GIPFELI_HEADERS:deps/gipfeli/include/%.h=phu/%.pcm)
PROTOBUF_HEADERS := $(wildcard deps/protobuf/src/google/protobuf/*.h) $(wildcard deps/protobuf/src/google/protobuf/*/*.h)
PROTOBUF_PRECOMPILED_HEADERS := $(PROTOBUF_HEADERS:deps/protobuf/src/%.h=phu/%.pcm)
GMOCK_PRECOMPILED_HEADERS := phu/gmock/gmock.pcm
GTEST_PRECOMPILED_HEADERS := phu/gtest/gtest.pcm

$(info GIPFELI_PRECOMPILED_HEADERS = $(GIPFELI_PRECOMPILED_HEADERS))

# Dependencies
pmu/principia/base/traits.pcm:	\
	pmu/std.pcm 
pmu/principia/base/not_null.pcm:                 \
	phu/absl/log/check.pcm			 \
	pmu/std.pcm                              \
	pmu/principia/base/not_constructible.pcm \
	pmu/principia/base/traits.pcm           

obj/base/not_null_test.o: \
	base/not_null_test.cpp \
	pmu/principia/base/not_null.pcm \
	phu/base/macros.pcm \
	phu/absl/hash/hash_testing.pcm \
	phu/gmock/gmock.pcm \
	phu/gtest/gtest.pcm \
	pmu/std.pcm \
	pmu/principia/base/not_null.pcm \
	pmu/principia/base/not_constructible.pcm \
	pmu/principia/base/traits.pcm \
	phu/absl/log/check.pcm  # Voidify::operator&& weirdness.

obj/base/not_null.o: \
	pmu/principia/base/not_null.pcm \
	phu/base/macros.pcm \
	phu/absl/hash/hash_testing.pcm \
	phu/gmock/gmock.pcm \
	phu/gtest/gtest.pcm \
	pmu/std.pcm \
	pmu/principia/base/not_null.pcm \
	pmu/principia/base/not_constructible.pcm \
	pmu/principia/base/traits.pcm 

obj/base/pull_serializer_test.o: \
	base/pull_serializer_test.cpp \
	pmu/std.pcm \
	phu/absl/log/check.pcm \
	phu/absl/log/log.pcm \
	phu/gmock/gmock.pcm \
	phu/gtest/gtest.pcm \
	phu/gipfeli/gipfeli.pcm \
	phu/gipfeli/compression.pcm \
	phu/serialization/geometry.pb.pcm \
	phu/serialization/quantities.pb.pcm \
	phu/serialization/physics.pb.pcm \
	pmu/principia/base/array.pcm \
	pmu/principia/base/not_null.pcm \
	pmu/principia/base/not_constructible.pcm \
	pmu/principia/base/traits.pcm \
	pmu/principia/base/sink_source.pcm \
	pmu/principia/base/pull_serializer.pcm

pmu/principia/base/pull_serializer.pcm: \
	phu/absl/log/check.pcm \
	phu/absl/base/thread_annotations.pcm \
	phu/absl/synchronization/mutex.pcm \
	phu/gipfeli/compression.pcm \
	phu/google/protobuf/io/zero_copy_stream.pcm \
	phu/google/protobuf/message.pcm \
	pmu/std.pcm \
	pmu/principia/base/not_null.pcm \
	pmu/principia/base/not_constructible.pcm \
	pmu/principia/base/traits.pcm \
	pmu/principia/base/sink_source.pcm \
	pmu/principia/base/array.pcm

pmu/principia/base/sink_source.pcm: \
	pmu/std.pcm \
	phu/gipfeli/sinksource.pcm \
	pmu/principia/base/array.pcm

pmu/principia/base/array.pcm: \
	phu/absl/log/check.pcm \
	phu/absl/log/log.pcm \
	pmu/std.pcm
bin/base/not_null_test: \
	obj/base/not_null_test.o\
	obj/gmock_main.o\
	obj/gmock-all.o\
	obj/gtest-all.o\
	obj/base/not_null.o

obj/base/array.o: \
	pmu/principia/base/array.pcm \
	pmu/std.pcm

obj/base/pull_serializer.o: \
	pmu/principia/base/pull_serializer.pcm \
	pmu/principia/base/not_constructible.pcm \
	pmu/principia/base/not_null.pcm \
	pmu/principia/base/traits.pcm \
	pmu/principia/base/array.pcm \
	pmu/principia/base/sink_source.pcm \
	pmu/std.pcm

obj/base/sink_source.o: \
	pmu/principia/base/sink_source.pcm \
	pmu/principia/base/array.pcm \
	pmu/std.pcm

bin/base/pull_serializer_test: \
	obj/base/pull_serializer_test.o \
	obj/base/array.o \
	obj/gmock_main.o \
	obj/gmock-all.o \
	obj/gtest-all.o \
	obj/base/pull_serializer.o \
	obj/base/sink_source.o \
	obj/base/not_null.o \
	$(SERIALIZATION_OBJECTS)

$(SERIALIZATION_OBJECTS): obj/%.o: %.cc

pmu/std.pcm: /usr/lib/llvm-23/share/libc++/v1/std.cppm
	@mkdir -p $(@D)
	clang++-23 $(ARGS) --precompile $< -o $@
	
phu/%.pcm:
	@mkdir -p $(@D)
	clang++-23 $(ARGS) -fmodule-header $< -o $@


# Note the -xc++-module because Clang does not know about the .ixx extension.
pmu/principia/%.pcm: %.ixx
	@mkdir -p $(@D)
	clang++-23 $(ARGS) \
		$(addprefix -fmodule-file=, $(filter phu/%.pcm, $^)) \
		$(foreach module, $(filter pmu/%.pcm, $^), -fmodule-file=$(subst /,.,$(module:pmu/%.pcm=%))=$(module)) \
		--precompile -xc++-module $< -o $@

$(PRINCIPIA_PRECOMPILED_HEADERS): phu/%.pcm: %.hpp
$(SERIALIZATION_PRECOMPILED_HEADERS): phu/%.pcm: %.h
$(ABSL_PRECOMPILED_HEADERS): phu/%.pcm: deps/abseil-cpp/install/include/%.h
$(GIPFELI_PRECOMPILED_HEADERS): phu/%.pcm: deps/gipfeli/include/%.h
$(GMOCK_PRECOMPILED_HEADERS): phu/%.pcm: deps/googletest/googlemock/include/%.h
$(GTEST_PRECOMPILED_HEADERS): phu/%.pcm: deps/googletest/googletest/include/%.h
$(PROTOBUF_PRECOMPILED_HEADERS): phu/%.pcm: deps/protobuf/src/%.h

obj/gmock_main.o: deps/googletest/googlemock/src/gmock_main.cc
	clang++-23 $(ARGS) --compile $< -o $@

obj/gmock-all.o: deps/googletest/googlemock/src/gmock-all.cc
	clang++-23 $(ARGS) --compile $< -o $@
obj/gtest-all.o: deps/googletest/googletest/src/gtest-all.cc
	clang++-23 $(ARGS) --compile $< -o $@

obj/%.o:
	@mkdir -p $(@D)
	clang++-23 $(ARGS) \
		$(addprefix -fmodule-file=, $(filter phu/%.pcm, $^)) \
		$(foreach module, $(filter pmu/%.pcm, $^), -fmodule-file=$(subst /,.,$(module:pmu/%.pcm=%))=$(module)) \
		--compile $< -o $@

bin/%: obj/%.o
	@mkdir -p $(@D)
	clang++-23 $(ARGS) $^ \
		deps/protobuf/libprotobuf.a \
		deps/protobuf/third_party/utf8_range/libutf8*.a \
		deps/gipfeli/libgipfeli.a \
		-Wl,--start-group deps/abseil-cpp/install/lib/libabsl_*.a -Wl,--end-group \
		deps/re2/libre2.a \
		-o $@
