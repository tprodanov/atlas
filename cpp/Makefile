CXX = g++
CXXFLAGS = -O3 -Wall -Werror -std=c++1y

EXE = atlas
BIN = bin
ARGS = argparse

OBJECTS = $(patsubst %.cpp,$(BIN)/%.o,$(wildcard *.cpp)) $(patsubst $(ARGS)/%.cpp,$(BIN)/$(ARGS)/%.o,$(wildcard $(ARGS)/*.cpp))

all: $(BIN)/$(ARGS) $(EXE)

$(EXE): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $(EXE)
	
$(BIN)/docopt/%.o: docopt/%.cpp
	$(CXX) $(CXXFLAGS) -c -MMD -o $@ $<
	
$(BIN)/%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -MMD -o $@ $<

$(BIN)/$(ARGS):
	mkdir -p $(BIN)/$(ARGS)

clean:
	rm -rf $(BIN)/*

.PHONY: clean all
