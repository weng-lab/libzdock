SHELL		= /bin/sh
PREFIX  = $(HOME)/local
CXX = g++-7
STRIP = strip
BIN_DIR = bin
BIN = $(BIN_DIR)/test
SRC = -I/home/vanderva/local/eigen3
BINOBJ = src/test.cpp
OBJ_DIR := build
LIB_DIR := lib
OPT		= -march=native -O3 -DEIGEN_USE_LAPACKE -Wall -pedantic
#OPT		= -march=native -O3 -DEIGEN_USE_LAPACKE -Wall -pedantic -funroll-loops
DEBUG		=
CXXFLAGS		= $(OPT) $(DEBUG)
LIBRARY		= pdb++
LIBARCH		= $(LIB_DIR)/lib$(LIBRARY).a
BIN_SOURCES = src/test.cpp
LIB_SOURCES = src/PDB.cpp src/pdb++.cpp src/pdb_chars.cpp src/pdb_read.cpp src/pdb_sprntf.cpp \
              src/pdb_sscanf.cpp src/pdb_type.cpp src/pdbinput.cpp src/ZDOCK.cpp src/Pruning.cpp \
              src/Constraints.cpp src/TransformLigand.cpp src/TransformMultimer.cpp src/TransformUtil.cpp

rwildcard=$(foreach d,$(wildcard $1*),$(call rwildcard,$d/,$2) $(filter $(subst *,%,$2),$d))
HEADERS = $(call rwildcard, src/, *.h) $(call rwildcard, src/, *.hpp) $(call rwildcard, src/, *.i)
OBJ := $(addprefix $(OBJ_DIR)/, $(patsubst %.cpp, %.o, $(call rwildcard, src/, *.cpp)))
BINOBJ = $(addprefix $(OBJ_DIR)/, $(patsubst %.cpp, %.o, $(BIN_SOURCES)))
LIBOBJ = $(addprefix $(OBJ_DIR)/, $(patsubst %.cpp, %.o, $(LIB_SOURCES)))

all:		Makefile $(OBJ_DIR) $(LIB_DIR) $(BIN_DIR) $(BIN)

$(LIBARCH): $(LIBOBJ)
	ar cru $(LIBARCH) $(LIBOBJ)

$(OBJ_DIR)/%.o: %.cpp $(HEADERS)
	@mkdir -p $(OBJ_DIR)/$(shell dirname $<)
	$(CXX) $(CXXFLAGS) $(SRC) -c -o $@ $<

$(BIN): $(BINOBJ) $(LIBARCH)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LD_FLAGS)
	$(STRIP) $@

$(OBJ_DIR) $(LIB_DIR) $(BIN_DIR):
	mkdir -p $@

clean:;		rm -f $(OBJ) $(LIBARCH) $(BIN)

