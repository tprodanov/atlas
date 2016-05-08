CXX = g++
CXXFLAGS = -O3 -Wall -Werror -std=c++1y

EXE = atlas
BIN = bin

OBJECTS = $(patsubst %.cpp,$(BIN)/%.o,$(wildcard *.cpp))

all: $(BIN) $(EXE)

$(EXE): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $(EXE)

$(BIN)/%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -MMD -o $@ $<

$(BIN):
	mkdir -p $(BIN)

clean:
	rm -rf $(BIN)/*

.PHONY: clean all
