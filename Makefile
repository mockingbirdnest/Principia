SOURCES=$(wildcard *.cpp) $(wildcard */*.cpp)
PROTO_SOURCES=$(wildcard */*.proto)

OBJECTS=$(SOURCES:.cpp=.o)
PROTO_HEADERS=$(PROTO_SOURCES:.proto=.pb.h)

BIN_DIR=bin
BIN=$(BIN_DIR)/principia

INCLUDE=-I. -Iglog/src -Iprotobuf/protobuf-2.6.1/src -Ibenchmark/include

CPPC=clang++
SHARED_ARGS=-std=c++1y -Ofast -g -ggdb -m64 -mmmx -msse -msse2 -m3dnow -fexceptions -ferror-limit=0 # -Wall -Wpedantic 
COMPILE_ARGS=-c $(SHARED_ARGS) $(INCLUDE)
LINK_ARGS=$(SHARED_ARGS)
LIBS=

$(BIN): $(PROTO_HEADERS) $(OBJECTS) Makefile $(BIN_DIR)
	$(CPPC) $(LINK_ARGS) $(OBJECTS) $(INCLUDE) $(LIBS) -o $(BIN)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

%.pb.h: %.proto
	protobuf/protobuf-2.6.1/src/protoc $< --cpp_out=.

%.o: %.cpp Makefile
	$(CPPC) $(COMPILE_ARGS) $< -MMD

run: $(BIN)
	./$(BIN)

clean:
	rm $(BIN) $(PROTO_HEADERS) $(OBJECTS); true
