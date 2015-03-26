SOURCES=$(wildcard *.cpp) $(wildcard quantities/*.cpp) $(wildcard base/*.cpp) $(wildcard integrators/*.cpp) $(wildcard ksp_plugin/*.cpp) $(wildcard geometry/*.cpp)
PROTO_SOURCES=$(wildcard */*.proto)

OBJECTS=$(SOURCES:.cpp=.o)
VERSIOON_HEADER=base/version.hpp
PROTO_HEADERS=$(PROTO_SOURCES:.proto=.pb.h)

BIN_DIR=bin
BIN=$(BIN_DIR)/principia.so

INCLUDE=-I. -I../glog/src -I../protobuf-2.6.1/src -I../benchmark/include -I../gmock-1.7.0/gtest/include -I../gmock-1.7.0/include

CPPC=clang++
SHARED_ARGS=-std=c++1y -stdlib=libc++ -O3 -g -ggdb -m64 -mmmx -msse -msse2 -m3dnow -fexceptions -ferror-limit=0 # -Wall -Wpedantic 
COMPILE_ARGS=-c $(SHARED_ARGS) $(INCLUDE)
LINK_ARGS=$(SHARED_ARGS)
LIBS=

$(BIN): $(VERSION_HEADER) $(PROTO_HEADERS) $(OBJECTS) Makefile $(BIN_DIR)
	$(CPPC) $(LINK_ARGS) $(OBJECTS) $(INCLUDE) $(LIBS) -o $(BIN)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(VERSION_HEADER): .git
	./generate_version_header.sh

%.pb.h: %.proto
	../protobuf/src/protoc $< --cpp_out=.

%.o: %.cpp Makefile
	$(CPPC) $(COMPILE_ARGS) $< -o $@ 

run: $(BIN)
	./$(BIN)

clean:
	rm $(BIN) $(PROTO_HEADERS) $(OBJECTS); true
