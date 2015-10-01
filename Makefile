CXX := clang++


VERSION_HEADER := base/version.hpp

CPP_SOURCES := ksp_plugin/plugin.cpp ksp_plugin/interface.cpp ksp_plugin/physics_bubble.cpp 
PROTO_SOURCES := $(wildcard */*.proto)
PROTO_CC_SOURCES := $(PROTO_SOURCES:.proto=.pb.cc)
PROTO_HEADERS := $(PROTO_SOURCES:.proto=.pb.h)

OBJECTS := $(CPP_SOURCES:.cpp=.o)
PROTO_OBJECTS := $(PROTO_CC_SOURCES:.cc=.o)
TEST_DIRS := base geometry integrators ksp_plugin_test physics quantities testing_utilities numerics
TEST_BINS := $(addsuffix /test,$(TEST_DIRS))

PROJECT_DIR := ksp_plugin_adapter/
SOLUTION_DIR := ./
ADAPTER_BUILD_DIR := ksp_plugin_adapter/obj
ADAPTER_CONFIGURATION := Debug
FINAL_PRODUCTS_DIR := Debug
ADAPTER := $(ADAPTER_BUILD_DIR)/$(ADAPTER_CONFIGURATION)/ksp_plugin_adapter.dll

LIB_DIR := $(FINAL_PRODUCTS_DIR)/GameData/Principia
LIB := $(LIB_DIR)/principia.so

DEP_DIR := deps
LIBS := $(DEP_DIR)/protobuf/src/.libs/libprotobuf.a $(DEP_DIR)/glog/.libs/libglog.a -lpthread -lc++ -lc++abi
TEST_INCLUDES := -I$(DEP_DIR)/googlemock/include -I$(DEP_DIR)/googletest/include -I $(DEP_DIR)/googlemock/ -I $(DEP_DIR)/googletest/ -I $(DEP_DIR)/eggsperimental_filesystem/
INCLUDES := -I. -I$(DEP_DIR)/glog/src -I$(DEP_DIR)/protobuf/src -I$(DEP_DIR)/benchmark/include -I$(DEP_DIR)/Optional $(TEST_INCLUDES)
SHARED_ARGS := -std=c++14 -stdlib=libc++ -O3 -g -fPIC -fexceptions -ferror-limit=0 -fno-omit-frame-pointer -Wall -Wpedantic \
	-DPROJECT_DIR='std::experimental::filesystem::path("$(PROJECT_DIR)")'\
	-DSOLUTION_DIR='std::experimental::filesystem::path("$(SOLUTION_DIR)")'

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

CXXFLAGS := -c $(SHARED_ARGS) $(INCLUDES)
LDFLAGS := $(SHARED_ARGS)


.PHONY: all adapter lib tests check plugin run_tests clean
.DEFAULT_GOAL := plugin

##### CONVENIENCE TARGETS #####
all: $(LIB) $(ADAPTER) tests

adapter: $(ADAPTER)
lib: $(LIB)

tests: $(TEST_BINS)

check: run_tests

##### CORE #####
$(ADAPTER):
	$(MDTOOL) build -c:$(ADAPTER_CONFIGURATION) ksp_plugin_adapter/ksp_plugin_adapter.csproj

$(LIB): $(VERSION_HEADER) $(PROTO_HEADERS) $(PROTO_OBJECTS) $(OBJECTS)
	$(CXX) -shared $(LDFLAGS) $(PROTO_OBJECTS) $(OBJECTS) -o $(LIB) $(LIBS) 

$(LIB_DIR):
	mkdir -p $(LIB_DIR)

$(VERSION_HEADER): .git
	./generate_version_header.sh

%.pb.cc %.pb.h: %.proto
	$(DEP_DIR)/protobuf/src/protoc $< --cpp_out=.

%.o: %.cpp 
	$(CXX) $(CXXFLAGS) $< -o $@ 

%.o: %.cc $(PROTO_HEADERS)
	$(CXX) $(CXXFLAGS) $< -o $@ 

##### DISTRIBUTION #####
plugin: $(ADAPTER) $(LIB)
	cd $(FINAL_PRODUCTS_DIR); zip -r Principia-$(UNAME_S)-$(shell git describe)-$(shell date "+%Y-%m-%d").zip GameData/

##### TESTS #####
run_tests: tests
	@echo "Cake, and grief counseling, will be available at the conclusion of the test."
	-base/test
	-geometry/test
	-integrators/test
	-ksp_plugin_test/test
	-physics/test
	-quantities/test
	-testing_utilities/test

TEST_LIBS=$(DEP_DIR)/protobuf/src/.libs/libprotobuf.a $(DEP_DIR)/glog/.libs/libglog.a -lpthread

GMOCK_SOURCE=$(DEP_DIR)/googlemock/src/gmock-all.cc $(DEP_DIR)/googlemock/src/gmock_main.cc $(DEP_DIR)/googletest/src/gtest-all.cc
GMOCK_OBJECTS=$(GMOCK_SOURCE:.cc=.o)

test_objects = $(patsubst %.cpp,%.o,$(wildcard $*/*.cpp))
ksp_plugin_test_objects = $(patsubst %.cpp,%.o,$(wildcard ksp_plugin/*.cpp)) $(patsubst %.cpp,%.o,$(wildcard ksp_plugin_test/*.cpp))

# We need to special-case ksp_plugin_test because it requires object files from ksp_plugin. The other tests don't do this.
.SECONDEXPANSION:
ksp_plugin_test/test: $$(ksp_plugin_test_objects) $(GMOCK_OBJECTS) $(PROTO_OBJECTS)
	$(CXX) $(LDFLAGS) $^ $(TEST_LIBS) -o $@

.SECONDEXPANSION:
%/test: $$(test_objects) $(GMOCK_OBJECTS) $(PROTO_OBJECTS)
	$(CXX) $(LDFLAGS) $^ $(TEST_LIBS) -o $@

##### CLEAN #####
clean_test-%:
	rm -f $(test_objects)

clean:  $(addprefix clean_test-,$(TEST_DIRS))
	rm -rf $(ADAPTER_BUILD_DIR) $(FINAL_PRODUCTS_DIR)
	rm -f $(LIB) $(VERSION_HEADER) $(PROTO_HEADERS) $(PROTO_CC_SOURCES) $(OBJECTS) $(PROTO_OBJECTS) $(TEST_BINS) $(LIB) $(ksp_plugin_test_objects)
