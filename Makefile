CXX := clang++


#TESTING CRAP REMOVE
b.test : a.test
	echo "b?"
	if test -e b.test; then echo "b exists"; else touch b.test; echo "made b"; fi
c.test : b.test
	echo "making c"
	touch c.test

PLUGIN_TRANSLATION_UNITS       := $(wildcard ksp_plugin/*.cpp)
PLUGIN_TEST_TRANSLATION_UNITS  := $(wildcard ksp_plugin_test/*.cpp)
JOURNAL_TRANSLATION_UNITS      := $(wildcard journal/*.cpp)
MOCK_TRANSLATION_UNITS         := $(wildcard */mock_*.cpp)
TEST_TRANSLATION_UNITS         := $(wildcard */*_test.cpp)
TEST_OR_MOCK_TRANSLATION_UNITS := $(TEST_TRANSLATION_UNITS) $(MOCK_TRANSLATION_UNITS)
TOOLS_TRANSLATION_UNITS        := $(wildcard tools/*.cpp)
LIBRARY_TRANSLATION_UNITS      := $(filter-out $(TEST_OR_MOCK_TRANSLATION_UNITS), $(wildcard */*.cpp))
JOURNAL_LIB_TRANSLATION_UNITS  := $(filter-out $(TEST_OR_MOCK_TRANSLATION_UNITS), $(wildcard journal/*.cpp))
BASE_LIB_TRANSLATION_UNITS     := $(filter-out $(TEST_OR_MOCK_TRANSLATION_UNITS), $(wildcard base/*.cpp))
PROTO_FILES                    := $(wildcard */*.proto)
PROTO_TRANSLATION_UNITS        := $(PROTO_FILES:.proto=.pb.cc)
PROTO_HEADERS                  := $(PROTO_FILES:.proto=.pb.h)

DEP_DIR := deps/

GMOCK_TRANSLATION_UNITS :=                     \
	$(DEP_DIR)googlemock/src/gmock-all.cc  \
	$(DEP_DIR)googlemock/src/gmock_main.cc \
	$(DEP_DIR)googletest/src/gtest-all.cc

VERSION_HEADER := base/version.generated.h

GENERATED_PROFILES :=                    \
	journal/profiles.generated.h     \
	journal/profiles.generated.cc    \
	journal/player.generated.cc      \
	ksp_plugin/interface.generated.h \
	ksp_plugin_adapter/interface.generated.cs

PROJECT_DIR := ksp_plugin_adapter/
SOLUTION_DIR := ./
ADAPTER_BUILD_DIR := ksp_plugin_adapter/obj
ADAPTER_CONFIGURATION := Debug
FINAL_PRODUCTS_DIR := Debug
ADAPTER := $(ADAPTER_BUILD_DIR)/$(ADAPTER_CONFIGURATION)/ksp_plugin_adapter.dll

LIB_DIR := $(FINAL_PRODUCTS_DIR)/GameData/Principia/Linux64
LIB := $(LIB_DIR)/principia.so
TEST_LIBS := $(DEP_DIR)benchmark/src/libbenchmark.a
LIBS      := $(DEP_DIR)/protobuf/src/.libs/libprotobuf.a $(DEP_DIR)/glog/.libs/libglog.a -lpthread -lc++ -lc++abi
TEST_INCLUDES := -I$(DEP_DIR)googlemock/include -I$(DEP_DIR)googletest/include -I$(DEP_DIR)googlemock/ -I$(DEP_DIR)googletest/ -I$(DEP_DIR)benchmark/include 
INCLUDES := -I. -I$(DEP_DIR)glog/src -I$(DEP_DIR)protobuf/src -I$(DEP_DIR)Optional -I$(DEP_DIR)eggsperimental_filesystem/
SHARED_ARGS := -std=c++14 -stdlib=libc++ -O3 -g -fPIC -fexceptions -ferror-limit=1 -fno-omit-frame-pointer -Wall -Wpedantic \
	-DPROJECT_DIR='std::experimental::filesystem::path("$(PROJECT_DIR)")'\
	-DSOLUTION_DIR='std::experimental::filesystem::path("$(SOLUTION_DIR)")' \
	-DNDEBUG

COMPILER_OPTIONS = -c $(SHARED_ARGS) $(INCLUDES)

# detect OS
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    UNAME_M := $(shell uname -m)
    ifeq ($(UNAME_M),x86_64)
        SHARED_ARGS += -m64
    else
        SHARED_ARGS += -m32
    endif
    MDTOOL := mdtool
endif
ifeq ($(UNAME_S),Darwin)
    SHARED_ARGS += -mmacosx-version-min=10.7 -arch i386
    MDTOOL ?= "/Applications/Xamarin Studio.app/Contents/MacOS/mdtool"
endif

########## Dependency resolution

BUILD_DIRECTORY := build/

TEST_OR_MOCK_DEPENDENCIES := $(addprefix $(BUILD_DIRECTORY), $(TEST_OR_MOCK_TRANSLATION_UNITS:.cpp=.d))
TOOLS_DEPENDENCIES        := $(addprefix $(BUILD_DIRECTORY), $(TOOLS_TRANSLATION_UNITS:.cpp=.d))
LIBRARY_DEPENDENCIES      := $(addprefix $(BUILD_DIRECTORY), $(LIBRARY_TEST_TRANSLATION_UNITS:.cpp=.d))
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
	sed 's!.*\.o[ :]*!$(BUILD_DIRECTORY)$*.o $@ : !g' < $@.temp > $@
	rm -f $@.temp

$(TEST_OR_MOCK_DEPENDENCIES): $(BUILD_DIRECTORY)%.d: %.cpp | $(PROTO_HEADERS)
	@mkdir -p $(@D)
	$(CXX) -M $(COMPILER_OPTIONS) $(TEST_INCLUDES) $< > $@.temp
	sed 's!.*\.o[ :]*!$(BUILD_DIRECTORY)$*.o $@ : !g' < $@.temp > $@
	rm -f $@.temp

include $(LIBRARY_DEPENDENCIES)
include $(TEST_DEPENDENCIES)

########## Compilation

##### C# and C++ code generation

$(VERSION_HEADER): .git
	./generate_version_header.sh

# We don't do dependency resolution on the protos; we compile them all at once.
$(PROTO_HEADERS) $(PROTO_TRANSLATION_UNITS): $(PROTO_FILES)
	$(DEP_DIR)/protobuf/src/protoc -I $(DEP_DIR)/protobuf/src/ -I . $^ --cpp_out=.

$(GENERATED_PROFILES) : $(TOOLS_BIN)
	$^ generate_profiles

##### C++ compilation

OBJ_DIRECTORY := obj/

TEST_OR_MOCK_OBJECTS := $(addprefix $(OBJ_DIRECTORY), $(TEST_OR_MOCK_TRANSLATION_UNITS:.cpp=.o))
LIBRARY_OBJECTS      := $(addprefix $(OBJ_DIRECTORY), $(LIBRARY_TRANSLATION_UNITS:.cpp=.o))
PROTO_OBJECTS        := $(addprefix $(OBJ_DIRECTORY), $(PROTO_TRANSLATION_UNITS:.cc=.o))
GMOCK_OBJECTS        := $(addprefix $(OBJ_DIRECTORY), $(GMOCK_TRANSLATION_UNITS:.cc=.o))
TOOLS_OBJECTS        := $(addprefix $(OBJ_DIRECTORY), $(TOOLS_TRANSLATION_UNITS:.cpp=.o))
PLUGIN_OBJECTS       := $(addprefix $(OBJ_DIRECTORY), $(PLUGIN_TRANSLATION_UNITS:.cpp=.o))
JOURNAL_LIB_OBJECTS  := $(addprefix $(OBJ_DIRECTORY), $(JOURNAL_LIB_TRANSLATION_UNITS:.cpp=.o))
BASE_LIB_OBJECTS     := $(addprefix $(OBJ_DIRECTORY), $(BASE_LIB_TRANSLATION_UNITS:.cpp=.o))
TEST_OBJECTS         := $(addprefix $(OBJ_DIRECTORY), $(TEST_TRANSLATION_UNITS:.cpp=.o))
MOCK_OBJECTS         := $(addprefix $(OBJ_DIRECTORY), $(MOCK_TRANSLATION_UNITS:.cpp=.o))

$(TEST_OR_MOCK_OBJECTS): $(OBJ_DIRECTORY)%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(COMPILER_OPTIONS) $(TEST_INCLUDES) $< -o $@

$(GMOCK_OBJECTS): $(OBJ_DIRECTORY)%.o: %.cc
	@mkdir -p $(@D)
	$(CXX) $(COMPILER_OPTIONS) $(TEST_INCLUDES) $< -o $@

$(LIBRARY_OBJECTS): $(OBJ_DIRECTORY)%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(COMPILER_OPTIONS) $< -o $@

$(PROTO_OBJECTS): $(OBJ_DIRECTORY)%.o: %.cc
	@mkdir -p $(@D)
	$(CXX) $(COMPILER_OPTIONS) $< -o $@

########## Linkage

BIN_DIRECTORY   := bin/

##### tools

TOOLS_BIN     := $(BIN_DIRECTORY)tools

$(TOOLS_BIN): $(TOOLS_OBJECTS) $(PROTO_OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(LDFLAGS) $^ $(LIBS) -o $@

##### KSP plugin

KSP_PLUGIN := $(BIN_DIRECTORY)principia.so

$(KSP_PLUGIN) : $(PROTO_OBJECTS) $(PLUGIN_OBJECTS) $(JOURNAL_LIB_OBJECTS) $(BASE_LIB_OBJECTS)
	@mkdir -p $(@D)
	$(CXX) -shared $(LDFLAGS) $^ $(LIBS) -o $@

CXXFLAGS := -c $(SHARED_ARGS) $(INCLUDES)
LDFLAGS := $(SHARED_ARGS)

##### Unit tests

TEST_BINS                    := $(addprefix $(BIN_DIRECTORY), $(TEST_TRANSLATION_UNITS:.cpp=))
PACKAGE_TEST_BINS            := $(addprefix $(BIN_DIRECTORY), $(addsuffix test, $(sort $(dir $(TEST_TRANSLATION_UNITS)))))
PLUGIN_DEPENDENT_TEST_BINS   := $(filter bin/ksp_plugin_test/% bin/journal/%, $(TEST_BINS))
PLUGIN_INDEPENDENT_TEST_BINS := $(filter-out $(PLUGIN_DEPENDENT_TEST_BINS), $(TEST_BINS))
PRINCIPIA_TEST_BIN           := $(BIN_DIRECTORY)test
log :
	@echo $(PACKAGE_TEST_BINS)

an_apple_pie : the_universe

$(PLUGIN_INDEPENDENT_TEST_BINS) : $(BIN_DIRECTORY)% : $(OBJ_DIRECTORY)%.o $(GMOCK_OBJECTS) $(PROTO_OBJECTS) $(BASE_LIB_OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(LDFLAGS) $^ $(LIBS) -o $@

# For tests that depend on the plugin, we link against the principia shared
# library instead of statically linking the objects.  Also note that we do not
# link the $(LIBS), since they are in the $(KSP_PLUGIN).  We still need pthread
# though.

$(PLUGIN_DEPENDENT_TEST_BINS) : $(BIN_DIRECTORY)% : $(OBJ_DIRECTORY)%.o $(MOCK_OBJECTS) $(GMOCK_OBJECTS) $(KSP_PLUGIN) $(TEST_LIBS)
	@mkdir -p $(@D)
	$(CXX) $(LDFLAGS) $^ -lpthread -o $@

$(PRINCIPIA_TEST_BIN) : $(TEST_OR_MOCK_OBJECTS) $(GMOCK_OBJECTS) $(KSP_PLUGIN) $(TEST_LIBS)
	@mkdir -p $(@D)
	$(CXX) $(LDFLAGS) $^ -lpthread -o $@

.PHONY: all adapter lib tests tools check plugin run_tests clean
.PRECIOUS: %.o $(PROTO_HEADERS) $(PROTO_CC_SOURCES) $(GENERATED_SOURCES)
.DEFAULT_GOAL := plugin
.SUFFIXES:

##### CONVENIENCE TARGETS #####
unit_tests : $(TEST_BINS)
all: $(LIB) $(ADAPTER) tests

adapter: $(ADAPTER)
lib: $(LIB)

tests: $(TEST_BINS)

tools: $(TOOLS_BIN)

check: run_tests

##### CORE #####
$(ADAPTER): $(GENERATED_SOURCES)
	$(MDTOOL) build -c:$(ADAPTER_CONFIGURATION) ksp_plugin_adapter/ksp_plugin_adapter.csproj

#$(TOOLS_BIN): $(PROTO_OBJECTS) $(TOOLS_OBJECTS) $(STATUS_OBJECTS)
#	$(CXX) $(LDFLAGS) $(PROTO_OBJECTS) $(TOOLS_OBJECTS) -o $(TOOLS_BIN) $(LIBS)

.SECONDEXPANSION:
$(LIB): $(PROTO_OBJECTS) $$(ksp_plugin_objects) $$(journal_objects) $(LIB_DIR) $(STATUS_OBJECTS)
	$(CXX) -shared $(LDFLAGS) $(PROTO_OBJECTS) $(STATUS_OBJECTS) $(ksp_plugin_objects) $(journal_objects) -o $(LIB) $(LIBS)

$(LIB_DIR):
	mkdir -p $(LIB_DIR)

$(GENERATED_SOURCES): $(TOOLS_BIN) serialization/journal.proto
	tools/tools generate_profiles

#%.pb.cc %.pb.h: %.proto
#	$(DEP_DIR)/protobuf/src/protoc -I $(DEP_DIR)/protobuf/src/ -I . $< --cpp_out=.

%.pb.o: %.pb.cc $(PROTO_HEADERS)
	$(CXX) $(CXXFLAGS) $< -o $@ 

tools/%.o: tools/%.cpp $(VERSION_HEADER) $(PROTO_HEADERS)
	$(CXX) $(CXXFLAGS) $< -o $@

%.o: %.cpp $(VERSION_HEADER) $(PROTO_HEADERS) $(GENERATED_SOURCES)
	$(CXX) $(CXXFLAGS) $< -o $@ 

%.o: %.cc $(VERSION_HEADER) $(PROTO_HEADERS) $(GENERATED_SOURCES)
	$(CXX) $(CXXFLAGS) $< -o $@ 

##### DISTRIBUTION #####
plugin: $(ADAPTER) $(LIB)
	cd $(FINAL_PRODUCTS_DIR); zip -r Principia-$(UNAME_S)-$(shell git describe)-$(shell date "+%Y-%m-%d").zip GameData/

##### TESTS #####
run_tests: tests
	@echo "Cake, and grief counseling, will be available at the conclusion of the test."
	-astronomy/test
	-base/test
	-geometry/test
	-integrators/test
	-ksp_plugin_test/test
	-physics/test
	-quantities/test
	-testing_utilities/test
	-numerics/test

#TEST_LIBS=$(DEP_DIR)/protobuf/src/.libs/libprotobuf.a $(DEP_DIR)/glog/.libs/libglog.a -lpthread

#GMOCK_SOURCE=$(DEP_DIR)/googlemock/src/gmock-all.cc $(DEP_DIR)/googlemock/src/gmock_main.cc $(DEP_DIR)/googletest/src/gtest-all.cc
#GMOCK_OBJECTS=$(GMOCK_SOURCE:.cc=.o)

test_objects = $(patsubst %.cpp,%.o,$(wildcard $(@D)/*.cpp))
#ksp_plugin_objects = $(patsubst %.cpp,%.o,$(wildcard ksp_plugin/*.cpp))
#journal_objects = journal/profiles.o journal/recorder.o

# We need to special-case ksp_plugin_test and journal because they require object files from ksp_plugin
# and journal.  The other tests don't do this.
#.SECONDEXPANSION:
#ksp_plugin_test/test: $$(ksp_plugin_objects) $$(journal_objects) $$(test_objects) $(GMOCK_OBJECTS) $(PROTO_OBJECTS) $(STATUS_OBJECTS)
#	$(CXX) $(LDFLAGS) $^ $(TEST_LIBS) -o $@

# We cannot link the player test because we do not have the benchmarks.  We only build the recorder test.
#.SECONDEXPANSION:
#journal/test: $$(ksp_plugin_objects) $$(journal_objects) journal/player.o journal/recorder_test.o $(GMOCK_OBJECTS) $(PROTO_OBJECTS) $(STATUS_OBJECTS)
#	$(CXX) $(LDFLAGS) $^ $(TEST_LIBS) -o $@

#.SECONDEXPANSION:
#%/test: $$(test_objects) $(GMOCK_OBJECTS) $(PROTO_OBJECTS) $(STATUS_OBJECTS)
#	$(CXX) $(LDFLAGS) $^ $(TEST_LIBS) -o $@

##### CLEAN #####
clean:
	rm -rf $(ADAPTER_BUILD_DIR) $(FINAL_PRODUCTS_DIR)
	rm -f $(LIB) $(VERSION_HEADER) $(PROTO_HEADERS) $(PROTO_CC_SOURCES) $(GENERATED_SOURCES) $(TEST_BINS) $(TOOLS_BIN) $(LIB) */*.o

##### EVERYTHING #####
# Compiles everything, but does not link anything.  Used to check standard compliance on code that we don't want to run on *nix.
compile_everything: $(patsubst %.cpp,%.o,$(wildcard */*.cpp))

##### IWYU #####
IWYU := deps/include-what-you-use/bin/include-what-you-use
IWYU_FLAGS := -Xiwyu --max_line_length=200 -Xiwyu --mapping_file="iwyu.imp" -Xiwyu --check_also=*/*.hpp
IWYU_NOSAFE_HEADERS := --nosafe_headers
REMOVE_BOM := for f in `ls */*.hpp && ls */*.cpp`; do awk 'NR==1{sub(/^\xef\xbb\xbf/,"")}1' $$f | awk NF{p=1}p > $$f.nobom; mv $$f.nobom $$f; done
RESTORE_BOM := for f in `ls */*.hpp && ls */*.cpp`; do awk 'NR==1{sub(/^/,"\xef\xbb\xbf\n")}1' $$f > $$f.withbom; mv $$f.withbom $$f; done
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

normalize_bom:
	$(REMOVE_BOM)
	$(RESTORE_BOM)
