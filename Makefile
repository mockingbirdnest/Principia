.SECONDEXPANSION:
PERCENT := %

# detect OS
UNAME_S := $(shell uname -s)
UNAME_M := $(shell uname -m)

CXX := clang++

VERSION_TRANSLATION_UNIT := base/version.generated.cc

PLUGIN_TRANSLATION_UNITS               := $(wildcard ksp_plugin/*.cpp)
PLUGIN_TEST_TRANSLATION_UNITS          := $(wildcard ksp_plugin_test/*.cpp)
JOURNAL_TRANSLATION_UNITS              := $(wildcard journal/*.cpp)
FAKE_OR_MOCK_TRANSLATION_UNITS         := $(wildcard */fake_*.cpp */mock_*.cpp)
BENCHMARK_TRANSLATION_UNITS            := $(wildcard benchmarks/*.cpp */benchmark.cpp)
TEST_TRANSLATION_UNITS                 := $(wildcard */*_test.cpp)
TEST_OR_FAKE_OR_MOCK_TRANSLATION_UNITS := $(TEST_TRANSLATION_UNITS) $(FAKE_OR_MOCK_TRANSLATION_UNITS)
TOOLS_TRANSLATION_UNITS                := $(wildcard tools/*.cpp)
LIBRARY_TRANSLATION_UNITS              := $(filter-out $(TEST_OR_FAKE_OR_MOCK_TRANSLATION_UNITS) $(BENCHMARK_TRANSLATION_UNITS), $(wildcard */*.cpp))
JOURNAL_LIB_TRANSLATION_UNITS          := $(filter-out $(TEST_OR_FAKE_OR_MOCK_TRANSLATION_UNITS), $(wildcard journal/*.cpp))
BASE_LIB_TRANSLATION_UNITS             := $(filter-out $(TEST_OR_FAKE_OR_MOCK_TRANSLATION_UNITS), $(wildcard base/*.cpp))
NUMERICS_LIB_TRANSLATION_UNITS         := $(filter-out $(TEST_OR_FAKE_OR_MOCK_TRANSLATION_UNITS), $(wildcard numerics/*.cpp))
PROTO_FILES                            := $(wildcard */*.proto)
PROTO_TRANSLATION_UNITS                := $(PROTO_FILES:.proto=.pb.cc)
PROTO_HEADERS                          := $(PROTO_FILES:.proto=.pb.h)

DEP_DIR := deps/

OBJ_DIRECTORY := obj/

BIN_DIRECTORY := bin/
TOOLS_BIN     := $(BIN_DIRECTORY)tools

GMOCK_TRANSLATION_UNITS := \
	$(DEP_DIR)googletest/googlemock/src/gmock-all.cc  \
	$(DEP_DIR)googletest/googletest/src/gtest-all.cc
GMOCK_MAIN_TRANSLATION_UNIT := $(DEP_DIR)googletest/googlemock/src/gmock_main.cc

GENERATED_PROFILES := \
	journal/profiles.generated.h     \
	journal/profiles.generated.cc    \
	journal/player.generated.cc      \
	ksp_plugin/interface.generated.h \
	ksp_plugin_adapter/interface.generated.cs

PROJECT_DIR           := ksp_plugin_adapter/
SOLUTION_DIR          := ./
ADAPTER_BUILD_DIR     := ksp_plugin_adapter/obj/
ADAPTER_CONFIGURATION := Release
FINAL_PRODUCTS_DIR    := Release/
ADAPTER               := $(ADAPTER_BUILD_DIR)$(ADAPTER_CONFIGURATION)/ksp_plugin_adapter.dll

ifeq ($(UNAME_S),Linux)
    PLUGIN_DIRECTORY      := $(FINAL_PRODUCTS_DIR)GameData/Principia/Linux64/
endif
ifeq ($(UNAME_S),Darwin)
    PLUGIN_DIRECTORY      := $(FINAL_PRODUCTS_DIR)GameData/Principia/MacOS64/
endif

TEST_LIBS     := $(DEP_DIR)benchmark/src/libbenchmark.a $(DEP_DIR)protobuf/src/.libs/libprotobuf.a
LIBS          := $(DEP_DIR)protobuf/src/.libs/libprotobuf.a \
	$(DEP_DIR)gipfeli/libgipfeli.a \
	$(DEP_DIR)abseil-cpp/absl/strings/libabsl_strings.a \
	$(DEP_DIR)abseil-cpp/absl/synchronization/libabsl_synchronization.a \
	$(DEP_DIR)abseil-cpp/absl/time/libabsl_*.a \
	$(DEP_DIR)abseil-cpp/absl/debugging/libabsl_*.a \
	$(DEP_DIR)abseil-cpp/absl/numeric/libabsl_*.a \
	$(DEP_DIR)abseil-cpp/absl/base/libabsl_*.a \
	$(DEP_DIR)glog/.libs/libglog.a -lpthread -lc++ -lc++abi
TEST_INCLUDES := \
	-I$(DEP_DIR)googletest/googlemock/include -I$(DEP_DIR)googletest/googletest/include \
	-I$(DEP_DIR)googletest/googlemock/ -I$(DEP_DIR)googletest/googletest/ -I$(DEP_DIR)benchmark/include
INCLUDES      := -I. -I$(DEP_DIR)glog/src -I$(DEP_DIR)protobuf/src -I$(DEP_DIR)compatibility/filesystem \
	-I$(DEP_DIR)gipfeli/include -I$(DEP_DIR)abseil-cpp
SHARED_ARGS   := \
	-std=c++1z -stdlib=libc++ -O3 -g                           \
	-fPIC -fexceptions -ferror-limit=1 -fno-omit-frame-pointer \
	-Wall -Wpedantic                                           \
	-DPROJECT_DIR='std::filesystem::path("$(PROJECT_DIR)")'    \
	-DSOLUTION_DIR='std::filesystem::path("$(SOLUTION_DIR)")'  \
	-DTEMP_DIR='std::filesystem::path("/tmp")'                 \
	-DNDEBUG

ifeq ($(UNAME_S),Linux)
    ifeq ($(UNAME_M),x86_64)
        SHARED_ARGS += -m64
    else
        SHARED_ARGS += -m32
    endif
    MDTOOL := mdtool
    LIBS += -lsupc++
    TEST_LIBS += -lsupc++
    SHAREDFLAG := -shared
endif
ifeq ($(UNAME_S),Darwin)
    INCLUDES += -I$(DEP_DIR)compatibility/optional -I$(DEP_DIR)Optional
    SHARED_ARGS += -mmacosx-version-min=10.12 -arch x86_64
    MDTOOL ?= "/Applications/Xamarin Studio.app/Contents/MacOS/mdtool"
    SHAREDFLAG := -dynamiclib
endif

COMPILER_OPTIONS := -c $(SHARED_ARGS) $(INCLUDES)
LDFLAGS := $(SHARED_ARGS)

########## Dependency resolution

BUILD_DIRECTORY := build/

TEST_OR_MOCK_DEPENDENCIES := $(addprefix $(BUILD_DIRECTORY), $(TEST_OR_FAKE_OR_MOCK_TRANSLATION_UNITS:.cpp=.d))
BENCHMARK_DEPENDENCIES    := $(addprefix $(BUILD_DIRECTORY), $(BENCHMARK_TRANSLATIONS_UNITS:.cpp=.d))
TOOLS_DEPENDENCIES        := $(addprefix $(BUILD_DIRECTORY), $(TOOLS_TRANSLATION_UNITS:.cpp=.d))
LIBRARY_DEPENDENCIES      := $(addprefix $(BUILD_DIRECTORY), $(LIBRARY_TRANSLATION_UNITS:.cpp=.d))
PLUGIN_DEPENDENCIES       := $(addprefix $(BUILD_DIRECTORY), $(PLUGIN_TRANSLATION_UNITS:.cpp=.d))
PLUGIN_TEST_DEPENDENCIES  := $(addprefix $(BUILD_DIRECTORY), $(PLUGIN_TEST_TRANSLATION_UNITS:.cpp=.d))
JOURNAL_DEPENDENCIES      := $(addprefix $(BUILD_DIRECTORY), $(JOURNAL_TRANSLATION_UNITS:.cpp=.d))

# As a prerequisite for listing the includes of things that depend on
# generated headers, we must generate said code.
# Note that the prerequisites for dependency files are order-only, for
# bootstrapping: once we have the dependency files, their own actual dependency
# on generated headers and on the translation unit will be listed there and
# recomputed as needed.
$(PLUGIN_DEPENDENCIES)              : | $(GENERATED_PROFILES)
$(PLUGIN_TEST_DEPENDENCIES)         : | $(GENERATED_PROFILES)
$(JOURNAL_DEPENDENCIES)             : | $(GENERATED_PROFILES)

$(LIBRARY_DEPENDENCIES): $(BUILD_DIRECTORY)%.d: %.cpp | $(PROTO_HEADERS)
	@mkdir -p $(@D)
	$(CXX) -M $(COMPILER_OPTIONS) $< > $@.temp
	sed 's!.*\.o[ :]*!$(OBJ_DIRECTORY)$*.o $@ : !g' < $@.temp > $@
	rm -f $@.temp

$(TEST_OR_MOCK_DEPENDENCIES): $(BUILD_DIRECTORY)%.d: %.cpp | $(PROTO_HEADERS)
	@mkdir -p $(@D)
	$(CXX) -M $(COMPILER_OPTIONS) $(TEST_INCLUDES) $< > $@.temp
	sed 's!.*\.o[ :]*!$(OBJ_DIRECTORY)$*.o $@ : !g' < $@.temp > $@
	rm -f $@.temp

ifneq ($(MAKECMDGOALS), clean)
include $(LIBRARY_DEPENDENCIES)
include $(TEST_OR_MOCK_DEPENDENCIES)
endif

########## Compilation

##### C# and C++ code generation

$(VERSION_TRANSLATION_UNIT): .git
	./generate_version_translation_unit.sh

# We don't do dependency resolution on the protos; we compile them all at once.
$(PROTO_HEADERS) $(PROTO_TRANSLATION_UNITS): $(PROTO_FILES)
	$(DEP_DIR)/protobuf/src/protoc -I $(DEP_DIR)/protobuf/src/ -I . $^ --cpp_out=.

$(GENERATED_PROFILES) : $(TOOLS_BIN)
	$^ generate_profiles

##### C++ compilation

TEST_OR_FAKE_OR_MOCK_OBJECTS := $(addprefix $(OBJ_DIRECTORY), $(TEST_OR_FAKE_OR_MOCK_TRANSLATION_UNITS:.cpp=.o))
LIBRARY_OBJECTS              := $(addprefix $(OBJ_DIRECTORY), $(LIBRARY_TRANSLATION_UNITS:.cpp=.o))
PROTO_OBJECTS                := $(addprefix $(OBJ_DIRECTORY), $(PROTO_TRANSLATION_UNITS:.cc=.o))
GMOCK_OBJECTS                := $(addprefix $(OBJ_DIRECTORY), $(GMOCK_TRANSLATION_UNITS:.cc=.o))
GMOCK_MAIN_OBJECT            := $(addprefix $(OBJ_DIRECTORY), $(GMOCK_MAIN_TRANSLATION_UNIT:.cc=.o))
BENCHMARK_OBJECTS            := $(addprefix $(OBJ_DIRECTORY), $(BENCHMARK_TRANSLATION_UNITS:.cpp=.o))
TOOLS_OBJECTS                := $(addprefix $(OBJ_DIRECTORY), $(TOOLS_TRANSLATION_UNITS:.cpp=.o))
PLUGIN_OBJECTS               := $(addprefix $(OBJ_DIRECTORY), $(PLUGIN_TRANSLATION_UNITS:.cpp=.o))
JOURNAL_LIB_OBJECTS          := $(addprefix $(OBJ_DIRECTORY), $(JOURNAL_LIB_TRANSLATION_UNITS:.cpp=.o))
VERSION_OBJECTS              := $(addprefix $(OBJ_DIRECTORY), $(VERSION_TRANSLATION_UNIT:.cc=.o))
BASE_LIB_OBJECTS             := $(addprefix $(OBJ_DIRECTORY), $(BASE_LIB_TRANSLATION_UNITS:.cpp=.o)) $(VERSION_OBJECTS)
NUMERICS_LIB_OBJECTS         := $(addprefix $(OBJ_DIRECTORY), $(NUMERICS_LIB_TRANSLATION_UNITS:.cpp=.o)) $(VERSION_OBJECTS)
TEST_OBJECTS                 := $(addprefix $(OBJ_DIRECTORY), $(TEST_TRANSLATION_UNITS:.cpp=.o))
FAKE_OR_MOCK_OBJECTS         := $(addprefix $(OBJ_DIRECTORY), $(FAKE_OR_MOCK_TRANSLATION_UNITS:.cpp=.o))

$(TEST_OR_FAKE_OR_MOCK_OBJECTS): $(OBJ_DIRECTORY)%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(COMPILER_OPTIONS) $(TEST_INCLUDES) $< -o $@

$(GMOCK_OBJECTS) $(GMOCK_MAIN_OBJECT): $(OBJ_DIRECTORY)%.o: %.cc
	@mkdir -p $(@D)
	$(CXX) $(COMPILER_OPTIONS) $(TEST_INCLUDES) $< -o $@

$(BENCHMARK_OBJECTS): $(OBJ_DIRECTORY)%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(COMPILER_OPTIONS) $(TEST_INCLUDES) $< -o $@

$(LIBRARY_OBJECTS): $(OBJ_DIRECTORY)%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(COMPILER_OPTIONS) $< -o $@

$(VERSION_OBJECTS): $(OBJ_DIRECTORY)%.o: %.cc
	@mkdir -p $(@D)
	$(CXX) $(COMPILER_OPTIONS) $< -o $@

$(PROTO_OBJECTS): $(OBJ_DIRECTORY)%.o: %.cc
	@mkdir -p $(@D)
	$(CXX) $(COMPILER_OPTIONS) $< -o $@

########## Linkage

##### tools

$(TOOLS_BIN): $(TOOLS_OBJECTS) $(PROTO_OBJECTS) $(NUMERICS_LIB_OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(LDFLAGS) $^ $(LIBS) -o $@

##### KSP plugin

KSP_PLUGIN := $(PLUGIN_DIRECTORY)principia.so

$(KSP_PLUGIN) : $(PROTO_OBJECTS) $(PLUGIN_OBJECTS) $(JOURNAL_LIB_OBJECTS) $(BASE_LIB_OBJECTS) $(NUMERICS_LIB_OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(SHAREDFLAG) $(LDFLAGS) $^ $(LIBS) -o $@

##### Tests

TEST_BINS                            := $(addprefix $(BIN_DIRECTORY), $(TEST_TRANSLATION_UNITS:.cpp=))
PACKAGE_TEST_BINS                    := $(addprefix $(BIN_DIRECTORY), $(addsuffix test, $(sort $(dir $(TEST_TRANSLATION_UNITS)))))
PLUGIN_DEPENDENT_TEST_BINS           := $(filter bin/ksp_plugin_test/% bin/journal/%, $(TEST_BINS))
PLUGIN_DEPENDENT_PACKAGE_TEST_BINS   := $(filter bin/ksp_plugin_test/% bin/journal/%, $(PACKAGE_TEST_BINS))
PLUGIN_INDEPENDENT_TEST_BINS         := $(filter-out $(PLUGIN_DEPENDENT_TEST_BINS), $(TEST_BINS))
PLUGIN_INDEPENDENT_PACKAGE_TEST_BINS := $(filter-out $(PLUGIN_DEPENDENT_PACKAGE_TEST_BINS), $(PACKAGE_TEST_BINS))
PRINCIPIA_TEST_BIN                   := $(BIN_DIRECTORY)test

$(TEST_BINS)          : $(BIN_DIRECTORY)% : $(OBJ_DIRECTORY)%.o
$(PACKAGE_TEST_BINS)  : $(BIN_DIRECTORY)%test : $$(filter $(OBJ_DIRECTORY)%$$(PERCENT), $(TEST_OBJECTS))
$(PRINCIPIA_TEST_BIN) : $(TEST_OBJECTS)

$(PLUGIN_INDEPENDENT_PACKAGE_TEST_BINS) $(PLUGIN_INDEPENDENT_TEST_BINS) : $(GMOCK_OBJECTS) $(GMOCK_MAIN_OBJECT) $(PROTO_OBJECTS) $(BASE_LIB_OBJECTS) $(NUMERICS_LIB_OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(LDFLAGS) $^ $(LIBS) -o $@

# For tests that depend on the plugin, we link against the principia shared
# library instead of statically linking the objects.
# NOTE(egg): this assumes that only the plugin-dependent tests need to be linked
# against mock objects.  The classes further up that are big enough to be mocked
# are likely to be highly templatized, so this will probably hold for a while.
$(PRINCIPIA_TEST_BIN) $(PLUGIN_DEPENDENT_PACKAGE_TEST_BINS) $(PLUGIN_DEPENDENT_TEST_BINS) : $(FAKE_OR_MOCK_OBJECTS) $(GMOCK_OBJECTS) $(GMOCK_MAIN_OBJECT) $(KSP_PLUGIN) $(BASE_LIB_OBJECTS) $(NUMERICS_LIB_OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(LDFLAGS) $^ $(TEST_LIBS) $(LIBS) -lpthread -o $@

########## Testing

TEST_TARGETS         := $(patsubst $(BIN_DIRECTORY)%, %, $(TEST_BINS))
PACKAGE_TEST_TARGETS := $(patsubst $(BIN_DIRECTORY)%, %, $(PACKAGE_TEST_BINS))

# make base/not_null_test compiles bin/base/not_null_test and runs it.
$(TEST_TARGETS) : % : $(BIN_DIRECTORY)%
	-$^

# make base/test compiles bin/base/test and runs it.
$(PACKAGE_TEST_TARGETS) : % : $(BIN_DIRECTORY)%
	-$^

test: $(PRINCIPIA_TEST_BIN)
	@echo "Cake, and grief counseling, will be available at the conclusion of the test."
	-$^

########## Benchmarks

PACKAGE_BENCHMARK_BINS := $(addprefix $(BIN_DIRECTORY), $(addsuffix benchmarks, $(sort $(dir $(BENCHMARK_TRANSLATION_UNITS)))))
PACKAGE_BENCHMARK_TARGET := $(patsubst $(BIN_DIRECTORY)%, %, $(PACKAGE_BENCHMARK_BINS))

PRINCIPIA_BENCHMARK_BIN := $(BIN_DIRECTORY)benchmark

$(PRINCIPIA_BENCHMARK_BIN) : $(BENCHMARK_OBJECTS) $(FAKE_OR_MOCK_OBJECTS) $(GMOCK_OBJECTS) $(KSP_PLUGIN) $(BASE_LIB_OBJECTS) $(NUMERICS_LIB_OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(LDFLAGS) $^ $(TEST_LIBS) -lpthread -o $@

benchmark: $(PRINCIPIA_BENCHMARK_BIN)
	-$^

########## Adapter

$(ADAPTER): $(GENERATED_PROFILES)
	$(MDTOOL) build -c:$(ADAPTER_CONFIGURATION) ksp_plugin_adapter/ksp_plugin_adapter.csproj

######### Distribution

release: $(ADAPTER) $(KSP_PLUGIN)
	cd $(FINAL_PRODUCTS_DIR); tar -c -z -f - GameData/ > principia_$(UNAME_S)-$(shell git describe --tags --always --dirty --abbrev=40 --long).tar.gz
########## Cleaning

clean:
	rm -rf $(BUILD_DIRECTORY) $(OBJ_DIRECTORY) $(BIN_DIRECTORY) $(ADAPTER_BUILD_DIR) $(FINAL_PRODUCTS_DIR)
	rm -f $(VERSION_TRANSLATION_UNIT) $(PROTO_TRANSLATION_UNITS) $(PROTO_HEADERS) $(GENERATED_PROFILES)

REMOVE_BOM := for f in `ls */*.hpp && ls */*.cpp`; do awk 'NR==1{sub(/^\xef\xbb\xbf/,"")}1' $$f | awk NF{p=1}p > $$f.nobom; mv $$f.nobom $$f; done
RESTORE_BOM := for f in `ls */*.hpp && ls */*.cpp`; do awk 'NR==1{sub(/^/,"\xef\xbb\xbf\n")}1' $$f > $$f.withbom; mv $$f.withbom $$f; done

normalize_bom:
	$(REMOVE_BOM)
	$(RESTORE_BOM)

##### clang-tidy

$(TEST_OR_FAKE_OR_MOCK_TRANSLATION_UNITS:.cpp=.cpp--tidy): %--tidy: %
	@mkdir -p $(@D)
	clang-tidy $< $(tidy_options) -- $(COMPILER_OPTIONS) $(TEST_INCLUDES)

$(LIBRARY_TRANSLATION_UNITS:.cpp=.cpp--tidy): %--tidy: %
	@mkdir -p $(@D)
	clang-tidy $< $(tidy_options) -- $(COMPILER_OPTIONS)

TIDY_TARGETS = $(TEST_OR_FAKE_OR_MOCK_TRANSLATION_UNITS:.cpp=.cpp--tidy) $(LIBRARY_TRANSLATION_UNITS:.cpp=.cpp--tidy)

########## Convenience targets
all: test release
tools: $(TOOLS_BIN)
adapter: $(ADAPTER)
plugin: $(KSP_PLUGIN)
each_test : $(TEST_TARGETS)
each_package_test : $(PACKAGE_TEST_TARGETS)
tidy : $(TIDY_TARGETS)

.PHONY: all tools adapter plugin each_test test release clean normalize_bom tidy $(TIDY_TARGETS) $(TEST_TARGETS) $(PACKAGE_TEST_TARGETS)
.PRECIOUS: %.o $(PROTO_HEADERS) $(PROTO_TRANSLATION_UNITS)
.DEFAULT_GOAL := all
.SUFFIXES:

an_apple_pie : the_universe

##### IWYU #####
IWYU := deps/include-what-you-use/bin/include-what-you-use
IWYU_FLAGS := -Xiwyu --max_line_length=200 -Xiwyu --mapping_file="iwyu.imp" -Xiwyu --check_also=*/*.hpp
IWYU_NOSAFE_HEADERS := --nosafe_headers
FIX_INCLUDES := deps/include-what-you-use/bin/fix_includes.py
IWYU_CHECK_ERROR := tee /dev/tty | test ! "`grep ' error: '`"
IWYU_TARGETS := $(wildcard */*.cpp)
IWYU_CLEAN := rm iwyu_generated_mappings.imp; rm */*.iwyu

iwyu_generate_mappings:
	{ ls */*_body.hpp && ls */*.generated.h; } | awk -f iwyu_generate_mappings.awk > iwyu_generated_mappings.imp

%.cpp!!iwyu: iwyu_generate_mappings
	$(IWYU) $(CXXFLAGS) $(subst !SLASH!,/, $*.cpp) $(IWYU_FLAGS) 2>&1 | tee $(subst !SLASH!,/, $*.iwyu) | $(IWYU_CHECK_ERROR)
	$(REMOVE_BOM)
	$(FIX_INCLUDES) < $(subst !SLASH!,/, $*.iwyu) | cat
	$(RESTORE_BOM)

iwyu: $(subst /,!SLASH!, $(addsuffix !!iwyu, $(IWYU_TARGETS)))
	$(IWYU_CLEAN)

%.cpp!!iwyu_unsafe: iwyu_generate_mappings
	$(IWYU) $(CXXFLAGS) $(subst !SLASH!,/, $*.cpp) $(IWYU_FLAGS) 2>&1 | tee $(subst !SLASH!,/, $*.iwyu) | $(IWYU_CHECK_ERROR)
	$(REMOVE_BOM)
	$(FIX_INCLUDES) $(IWYU_NOSAFE_HEADERS) < $(subst !SLASH!,/, $*.iwyu) | cat
	$(RESTORE_BOM)

iwyu_unsafe: $(subst /,!SLASH!, $(addsuffix !!iwyu_unsafe, $(IWYU_TARGETS)))
	$(IWYU_CLEAN)

iwyu_clean:
	$(IWYU_CLEAN)
