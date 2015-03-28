CPP_SOURCES=ksp_plugin/plugin.cpp ksp_plugin/interface.cpp ksp_plugin/physics_bubble.cpp 
#CPP_SOURCES=$(wildcard testing_utilities/*.cpp) 
#CPP_SOURCES=physics/n_body_system_test.cpp
#$(wildcard quantities/*.cpp) $(wildcard base/*.cpp) $(wildcard integrators/*.cpp) $(wildcard ksp_plugin/*.cpp) $(wildcard geometry/*.cpp) $(wildcard physics/*.cpp) $(wildcard benchmarks/*.cpp)
PROTO_SOURCES=$(wildcard */*.proto)
PROTO_CC_SOURCES=$(wildcard serialization/*.cc)
PROTO_OBJECTS=$(PROTO_CC_SOURCES:.cc=.o)

OBJECTS=$(CPP_SOURCES:.cpp=.o)
VERSION_HEADER=base/version.hpp
PROTO_HEADERS=$(PROTO_SOURCES:.proto=.pb.h)

LIB_DIR=lib
LIB=$(LIB_DIR)/principia.so

DEP_DIR=deps

INCLUDE=-I. -I$(DEP_DIR)/glog/src -I$(DEP_DIR)/protobuf/src -I$(DEP_DIR)/benchmark/include -I$(DEP_DIR)/gmock-1.7.0/gtest/include -I$(DEP_DIR)/gmock-1.7.0/include

CPPC=clang++
SHARED_ARGS=-std=c++1y -stdlib=libc++ -O3 -g -ggdb -m64 -fPIC -fexceptions -ferror-limit=0 # -Wall -Wpedantic 
COMPILE_ARGS=-c $(SHARED_ARGS) $(INCLUDE)
LINK_ARGS=$(SHARED_ARGS) 
LIB_PATHS=-L$(DEP_DIR)/glog/.libs/ -L$(DEP_DIR)/benchmark/src/ -L$(DEP_DIR)/protobuf/src/.libs/ -L$(DEP_DIR)/gmock-1.7.0/lib/.libs/
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
run_tests: base/test geometry/test integrators/test ksp_plugin_test/test physics/test quantities/test testing_utilities/test
	base/test
	geometry/test
	integrators/test
	ksp_plugin_test/test
	physics/test
	quantities/test
	testing_utilities/test

GMOCK_SOURCE=$(wildcard $(DEP_DIR)/gmock-1.7.0/fused-src/*.cc)
GMOCK_OBJECTS=$(GMOCK_SOURCE:.cc=.o)

BASE_TEST_SOURCES=$(wildcard base/*.cpp)
BASE_TEST_OBJECTS=$(BASE_TEST_SOURCES:.cpp=.o)
base/test: $(GMOCK_OBJECTS) $(BASE_TEST_OBJECTS) $(PROTO_OBJECTS) Makefile
	$(CPPC) $(LINK_ARGS) $(BASE_TEST_OBJECTS) $(GMOCK_OBJECTS) $(PROTO_OBJECTS) $(INCLUDE) $(LIB_PATHS) $(LIBS) -o $@

GEOMETRY_TEST_SOURCES=$(wildcard geometry/*.cpp)
GEOMETRY_TEST_OBJECTS=$(GEOMETRY_TEST_SOURCES:.cpp=.o)
geometry/test: $(GMOCK_OBJECTS) $(GEOMETRY_TEST_OBJECTS) $(PROTO_OBJECTS) Makefile
	$(CPPC) $(LINK_ARGS) $(GEOMETRY_TEST_OBJECTS) $(GMOCK_OBJECTS) $(PROTO_OBJECTS) $(INCLUDE) $(LIB_PATHS) $(LIBS) -o $@

INTEGRATOR_TEST_SOURCES=$(wildcard integrators/*.cpp)
INTEGRATOR_TEST_OBJECTS=$(INTEGRATOR_TEST_SOURCES:.cpp=.o)
integrators/test: $(GMOCK_OBJECTS) $(INTEGRATOR_TEST_OBJECTS) $(PROTO_OBJECTS) Makefile
	$(CPPC) $(LINK_ARGS) $(INTEGRATOR_TEST_OBJECTS) $(GMOCK_OBJECTS) $(PROTO_OBJECTS) $(INCLUDE) $(LIB_PATHS) $(LIBS) -o $@

PLUGIN_TEST_SOURCES=$(wildcard ksp_plugin_test/*.cpp)
PLUGIN_TEST_OBJECTS=$(PLUGIN_TEST_SOURCES:.cpp=.o)
ksp_plugin_test/test: $(GMOCK_OBJECTS) $(PLUGIN_TEST_OBJECTS) $(PROTO_OBJECTS) Makefile
	$(CPPC) $(LINK_ARGS) $(PLUGIN_TEST_OBJECTS) $(GMOCK_OBJECTS) $(PROTO_OBJECTS) $(INCLUDE) $(LIB_PATHS) $(LIBS) -o $@

PHYSICS_TEST_SOURCES=$(wildcard physics/*.cpp)
PHYSICS_TEST_OBJECTS=$(PHYSICS_TEST_SOURCES:.cpp=.o)
physics/test: $(GMOCK_OBJECTS) $(PHYSICS_TEST_OBJECTS) $(PROTO_OBJECTS) Makefile
	$(CPPC) $(LINK_ARGS) $(PHYSICS_TEST_OBJECTS) $(GMOCK_OBJECTS) $(PROTO_OBJECTS) $(INCLUDE) $(LIB_PATHS) $(LIBS) -o $@

QUANTITIES_TEST_SOURCES=$(wildcard quantities/*.cpp)
QUANTITIES_TEST_OBJECTS=$(QUANTITIES_TEST_SOURCES:.cpp=.o)
quantities/test: $(GMOCK_OBJECTS) $(QUANTITIES_TEST_OBJECTS) $(PROTO_OBJECTS) Makefile
	$(CPPC) $(LINK_ARGS) $(QUANTITIES_TEST_OBJECTS) $(GMOCK_OBJECTS) $(PROTO_OBJECTS) $(INCLUDE) $(LIB_PATHS) $(LIBS) -o $@

TESTING_UTILITIES_TEST_SOURCES=$(wildcard testing_utilities/*.cpp)
TESTING_UTILITIES_TEST_OBJECTS=$(TESTING_UTILITIES_TEST_SOURCES:.cpp=.o)
testing_utilities/test: $(GMOCK_OBJECTS) $(TESTING_UTILITIES_TEST_OBJECTS) $(PROTO_OBJECTS) Makefile
	$(CPPC) $(LINK_ARGS) $(TESTING_UTILITIES_TEST_OBJECTS) $(GMOCK_OBJECTS) $(PROTO_OBJECTS) $(INCLUDE) $(LIB_PATHS) $(LIBS) -o $@

clean:
	rm $(LIB) $(PROTO_HEADERS) $(OBJECTS) $(PROTO_OBJECTS) $(TEST_OBJS); true
