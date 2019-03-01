SHELL		= /bin/sh
PREFIX  = $(HOME)/local
CXX = clang++
STRIP = strip
BIN_DIR = bin
OBJ_DIR := build
LIB_DIR := lib
SRC = -Icontrib/eigen
OPT		= -march=native -O3 -DEIGEN_USE_LAPACKE -Wall -pedantic
#OPT		= -march=native -O3 -DEIGEN_USE_LAPACKE -Wall -pedantic -funroll-loops
DEBUG		=
CXXFLAGS		= $(OPT) $(DEBUG)

BINS = $(BIN_DIR)/createlig $(BIN_DIR)/createmultimer $(BIN_DIR)/pruning \
       $(BIN_DIR)/constraints $(BIN_DIR)/centroids $(BIN_DIR)/zdsplit \
       $(BIN_DIR)/zdunsplit
LIBRARY		= zdock
LIBARCH		= $(LIB_DIR)/lib$(LIBRARY).a
PYTHON_DIR = python
PYTHON_CLEAN = $(PYTHON_DIR)/build $(PYTHON_DIR)/dist $(PYTHON_DIR)/zdock.egg-info

LIB_SOURCES = src/libpdb++/pdbinput.cpp src/libpdb++/pdb_read.cpp src/libpdb++/pdb++.cpp \
             src/libpdb++/pdb_sscanf.cpp src/libpdb++/pdb_type.cpp src/libpdb++/pdb_sprntf.cpp \
             src/libpdb++/pdb_chars.cpp src/zdock/TransformMultimer.cpp \
             src/zdock/Constraints.cpp src/zdock/TransformLigand.cpp src/zdock/TransformUtil.cpp \
             src/zdock/ZDOCK.cpp src/pdb/PDB.cpp
INCLUDE_PATHS = -Isrc/libpdb++ -Isrc/zdock -Isrc/common -Isrc/pdb -Iinclude
SRC += $(INCLUDE_PATHS)
rwildcard=$(foreach d,$(wildcard $1*),$(call rwildcard,$d/,$2) $(filter $(subst *,%,$2),$d))
HEADERS = $(call rwildcard, include/, *.hpp) \
          $(call rwildcard, include/, *.h) \
          $(call rwildcard, src/, */*.hpp) \
          $(call rwildcard, src/, *.hpp) \
          $(call rwildcard, src/libpdb++, *.i)
OBJ := $(addprefix $(OBJ_DIR)/, $(patsubst %.cpp, %.o, $(call rwildcard, src/, *.cpp)))
LIBOBJ = $(addprefix $(OBJ_DIR)/, $(patsubst %.cpp, %.o, $(LIB_SOURCES)))

all:		Makefile $(OBJ_DIR) $(LIB_DIR) $(BIN_DIR) $(BINS)

$(LIBARCH): $(LIBOBJ)
	ar cru $(LIBARCH) $(LIBOBJ)

$(OBJ_DIR)/%.o: %.cpp $(HEADERS)
	@mkdir -p $(OBJ_DIR)/$(shell dirname $<)
	$(CXX) $(CXXFLAGS) $(SRC) -c -o $@ $<

$(BIN_DIR)/createlig: build/src/CreateLigand.o $(LIBARCH)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LD_FLAGS)
	$(STRIP) $@

$(BIN_DIR)/createmultimer: build/src/CreateMultimer.o $(LIBARCH)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LD_FLAGS)
	$(STRIP) $@

$(BIN_DIR)/pruning: build/src/Pruning.o $(LIBARCH)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LD_FLAGS)
	$(STRIP) $@

$(BIN_DIR)/constraints: build/src/FilterConstraints.o $(LIBARCH)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LD_FLAGS)
	$(STRIP) $@

$(BIN_DIR)/centroids: build/src/Centroids.o $(LIBARCH)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LD_FLAGS)
	$(STRIP) $@

$(BIN_DIR)/zdsplit: build/src/Split.o $(LIBARCH)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LD_FLAGS)
	$(STRIP) $@

$(BIN_DIR)/zdunsplit: build/src/UnSplit.o $(LIBARCH)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LD_FLAGS)
	$(STRIP) $@

$(OBJ_DIR) $(LIB_DIR) $(BIN_DIR):
	mkdir -p $@

clean:;		rm -Rf $(OBJ) $(LIBARCH) $(BINS) $(PYTHON_CLEAN)

