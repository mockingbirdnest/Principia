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

TEST_BIN=principia_tests
TEST_SOURCE=$(wildcard base/*test.cpp) $(wildcard geometry/*test.cpp) $(wildcard integrators/*test.cpp) $(wildcard quantities/*test.cpp)
TEST_OBJS=$(TEST_SOURCE:.cpp=.o)
GMOCK_SOURCE=$(wildcard $(DEP_DIR)/gmock-1.7.0/fused-src/*.cc)
GMOCK_OBJS=$(GMOCK_SOURCE:.cc=.o)

$(TEST_BIN): $(GMOCK_OBJS) $(TEST_OBJS) $(PROTO_OBJECTS) Makefile
	$(CPPC) $(LINK_ARGS)$(TEST_OBJS) $(GMOCK_OBJS) $(PROTO_OBJECTS) $(INCLUDE) $(LIB_PATHS) $(LIBS) -o $(TEST_BIN)

run_tests: $(TEST_BIN)
	./$(TEST_BIN)

clean:
	rm $(LIB) $(PROTO_HEADERS) $(OBJECTS) $(PROTO_OBJECTS) $(TEST_OBJS); true
