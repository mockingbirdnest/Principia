CPP_SOURCES=ksp_plugin/plugin.cpp ksp_plugin/interface.cpp ksp_plugin/physics_bubble.cpp 
PROTO_SOURCES=$(wildcard */*.proto)
PROTO_CC_SOURCES=$(wildcard serialization/*.cc)
PROTO_OBJECTS=$(PROTO_CC_SOURCES:.cc=.o)

OBJECTS=$(CPP_SOURCES:.cpp=.o)
VERSION_HEADER=base/version.hpp
PROTO_HEADERS=$(PROTO_SOURCES:.proto=.pb.h)

LIB_DIR=lib
LIB=$(LIB_DIR)/principia.so

DEP_DIR=deps

TEST_INCLUDE=-I$(DEP_DIR)/gmock/include -I$(DEP_DIR)/gtest/include -I$(DEP_DIR)/gmock -I$(DEP_DIR)/gtest
INCLUDE=-I. -I$(DEP_DIR)/glog/src -I$(DEP_DIR)/protobuf/src -I$(DEP_DIR)/benchmark/include $(TEST_INCLUDE)

CPPC=clang++
SHARED_ARGS=-std=c++1y -stdlib=libc++ -O3 -g -m64 -fPIC -fexceptions -ferror-limit=0 -fno-omit-frame-pointer # -Wall -Wpedantic 
COMPILE_ARGS=-c $(SHARED_ARGS) $(INCLUDE)
LINK_ARGS=$(SHARED_ARGS) 
LIB_PATHS=-L$(DEP_DIR)/glog/.libs/ -L$(DEP_DIR)/benchmark/src/ -L$(DEP_DIR)/protobuf/src/.libs/ -L$(DEP_DIR)/gmock/lib/.libs/
LIBS=-l:libc++.a -l:libprotobuf.a -l:libglog.a -lpthread

all: $(DEP_DIR) $(LIB) run_tests

$(LIB): $(VERSION_HEADER) $(PROTO_HEADERS) $(OBJECTS) $(PROTO_OBJECTS) Makefile $(LIB_DIR)
	$(CPPC) -shared $(LINK_ARGS) $(OBJECTS) $(PROTO_OBJECTS) $(INCLUDE) -o $(LIB) $(LIB_PATHS) $(LIBS) 

$(LIB_DIR):
	mkdir -p $(LIB_DIR)

$(DEP_DIR):
	./install_deps.sh

$(VERSION_HEADER): .git
	./generate_version_header.sh

%.pb.h: %.proto
	$(DEP_DIR)/protobuf/src/protoc $< --cpp_out=.

%.o: %.cpp 
	$(CPPC) $(COMPILE_ARGS) $< -o $@ 

%.o: %.cc $(PROTO_HEADERS)
	$(CPPC) $(COMPILE_ARGS) $< -o $@ 

##### TESTS #####
TEST_BINS=base/test geometry/test integrators/test ksp_plugin_test/test physics/test quantities/test testing_utilities/test
run_tests: $(TEST_BINS)
	base/test; true
	geometry/test; true
	integrators/test; true
	ksp_plugin_test/test; true
	physics/test; true
	quantities/test; true
	testing_utilities/test; true

TEST_LIBS=-l:libc++.a -l:libprotobuf.a -l:libglog.a -lpthread

GMOCK_SOURCE=$(DEP_DIR)/gmock/src/gmock-all.cc $(DEP_DIR)/gmock/src/gmock_main.cc $(DEP_DIR)/gtest/src/gtest-all.cc
GMOCK_OBJECTS=$(GMOCK_SOURCE:.cc=.o)

BASE_TEST_SOURCES=$(wildcard base/*.cpp)
BASE_TEST_OBJECTS=$(BASE_TEST_SOURCES:.cpp=.o)
base/test: $(GMOCK_OBJECTS) $(BASE_TEST_OBJECTS) $(PROTO_OBJECTS) Makefile
	$(CPPC) $(LINK_ARGS) $(BASE_TEST_OBJECTS) $(GMOCK_OBJECTS) $(PROTO_OBJECTS) $(INCLUDE) $(LIB_PATHS) $(TEST_LIBS) -o $@

GEOMETRY_TEST_SOURCES=$(wildcard geometry/*.cpp)
GEOMETRY_TEST_OBJECTS=$(GEOMETRY_TEST_SOURCES:.cpp=.o)
geometry/test: $(GMOCK_OBJECTS) $(GEOMETRY_TEST_OBJECTS) $(PROTO_OBJECTS) Makefile
	$(CPPC) $(LINK_ARGS) $(GEOMETRY_TEST_OBJECTS) $(GMOCK_OBJECTS) $(PROTO_OBJECTS) $(INCLUDE) $(LIB_PATHS) $(TEST_LIBS) -o $@

INTEGRATOR_TEST_SOURCES=$(wildcard integrators/*.cpp)
INTEGRATOR_TEST_OBJECTS=$(INTEGRATOR_TEST_SOURCES:.cpp=.o)
integrators/test: $(GMOCK_OBJECTS) $(INTEGRATOR_TEST_OBJECTS) $(PROTO_OBJECTS) Makefile
	$(CPPC) $(LINK_ARGS) $(INTEGRATOR_TEST_OBJECTS) $(GMOCK_OBJECTS) $(PROTO_OBJECTS) $(INCLUDE) $(LIB_PATHS) $(TEST_LIBS) -o $@

PLUGIN_TEST_SOURCES=$(wildcard ksp_plugin_test/*.cpp) $(wildcard ksp_plugin/*.cpp)
PLUGIN_TEST_OBJECTS=$(PLUGIN_TEST_SOURCES:.cpp=.o)
ksp_plugin_test/test: $(GMOCK_OBJECTS) $(PLUGIN_TEST_OBJECTS) $(PROTO_OBJECTS) Makefile
	$(CPPC) $(LINK_ARGS) $(PLUGIN_TEST_OBJECTS) $(GMOCK_OBJECTS) $(PROTO_OBJECTS) $(INCLUDE) $(LIB_PATHS) $(TEST_LIBS) -o $@

PHYSICS_TEST_SOURCES=$(wildcard physics/*.cpp)
PHYSICS_TEST_OBJECTS=$(PHYSICS_TEST_SOURCES:.cpp=.o)
physics/test: $(GMOCK_OBJECTS) $(PHYSICS_TEST_OBJECTS) $(PROTO_OBJECTS) Makefile
	$(CPPC) $(LINK_ARGS) $(PHYSICS_TEST_OBJECTS) $(GMOCK_OBJECTS) $(PROTO_OBJECTS) $(INCLUDE) $(LIB_PATHS) $(TEST_LIBS) -o $@

QUANTITIES_TEST_SOURCES=$(wildcard quantities/*.cpp)
QUANTITIES_TEST_OBJECTS=$(QUANTITIES_TEST_SOURCES:.cpp=.o)
quantities/test: $(GMOCK_OBJECTS) $(QUANTITIES_TEST_OBJECTS) $(PROTO_OBJECTS) Makefile
	$(CPPC) $(LINK_ARGS) $(QUANTITIES_TEST_OBJECTS) $(GMOCK_OBJECTS) $(PROTO_OBJECTS) $(INCLUDE) $(LIB_PATHS) $(TEST_LIBS) -o $@

TESTING_UTILITIES_TEST_SOURCES=$(wildcard testing_utilities/*.cpp)
TESTING_UTILITIES_TEST_OBJECTS=$(TESTING_UTILITIES_TEST_SOURCES:.cpp=.o)
testing_utilities/test: $(GMOCK_OBJECTS) $(TESTING_UTILITIES_TEST_OBJECTS) $(PROTO_OBJECTS) Makefile
	$(CPPC) $(LINK_ARGS) $(TESTING_UTILITIES_TEST_OBJECTS) $(GMOCK_OBJECTS) $(PROTO_OBJECTS) $(INCLUDE) $(LIB_PATHS) $(TEST_LIBS) -o $@

clean:
	rm $(LIB) $(PROTO_HEADERS) $(OBJECTS) $(PROTO_OBJECTS) $(BASE_TEST_OBJECTS) $(GEOMETRY_TEST_OBJECTS) $(INTEGRATOR_TEST_OBJECTS) $(PLUGIN_TEST_OBJECTS) $(PHYSICS_TEST_OBJS) $(QUANTITIES_TEST_OBJECTS) $(TESTING_UTILITIES_TEST_OBJECTS) $(TEST_BINS); true
