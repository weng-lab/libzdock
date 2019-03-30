SHELL		= /bin/sh
PREFIX  = $(HOME)/local
DOXYGEN = doxygen
CXX = g++-7
STRIP = strip
BIN_DIR = bin
OBJ_DIR := build
LIB_DIR := lib
DOC_DIR := doc
TEST_DIR := test
SRC = -Icontrib/eigen
OPT		= -march=native -O3 -DEIGEN_USE_LAPACKE -Wall -pedantic
TEST_SRC = -Icontrib/Catch2/single_include
TEST_OPT = -DDATADIR=$(realpath $(TEST_DIR)/data)
DEBUG		=
override CXXFLAGS	+= $(OPT) $(DEBUG)

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
TEST_SOURCES = $(call rwildcard, test/, *.cpp)
INCLUDE_PATHS = -Isrc/libpdb++ -Isrc/zdock -Isrc/common -Isrc/pdb -Iinclude
SRC += $(INCLUDE_PATHS)
rwildcard=$(foreach d,$(wildcard $1*),$(call rwildcard,$d/,$2) $(filter $(subst *,%,$2),$d))
HEADERS = $(call rwildcard, include/, *.hpp) \
          $(call rwildcard, include/, *.h) \
          $(call rwildcard, src/, */*.hpp) \
          $(call rwildcard, src/, *.hpp) \
          $(call rwildcard, test/, *.hpp) \
          $(call rwildcard, src/libpdb++, *.i)
OBJ := $(addprefix $(OBJ_DIR)/, $(patsubst %.cpp, %.o, $(call rwildcard, src/, *.cpp))) \
       $(addprefix $(OBJ_DIR)/, $(patsubst %.cpp, %.o, $(call rwildcard, test/, *.cpp)))
LIBOBJ = $(addprefix $(OBJ_DIR)/, $(patsubst %.cpp, %.o, $(LIB_SOURCES)))
TESTOBJ = $(addprefix $(OBJ_DIR)/, $(patsubst %.cpp, %.o, $(TEST_SOURCES)))

all:		Makefile $(OBJ_DIR) $(LIB_DIR) $(BIN_DIR) $(BINS)

test: Makefile $(TEST_DIR)/test

$(LIBARCH): $(LIBOBJ)
	ar cru $(LIBARCH) $(LIBOBJ)

$(OBJ_DIR)/src/%.o: src/%.cpp $(HEADERS)
	@mkdir -p $(OBJ_DIR)/$(shell dirname $<)
	$(CXX) $(CXXFLAGS) $(SRC) -c -o $@ $<

$(OBJ_DIR)/test/%.o: test/%.cpp $(HEADERS)
	@mkdir -p $(OBJ_DIR)/$(shell dirname $<)
	$(CXX) $(CXXFLAGS) $(TEST_OPT) $(TEST_SRC) $(SRC) -c -o $@ $<

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

$(TEST_DIR)/test: $(TESTOBJ) $(LIBARCH)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LD_FLAGS)
	$(STRIP) $@
	$(TEST_DIR)/test

$(OBJ_DIR) $(LIB_DIR) $(BIN_DIR):
	mkdir -p $@

doc/doxygen/html/index.html: Makefile doc/Doxyfile $(HEADERS) $(call rwildcard, src/, *.cpp) $(call rwildcard, test/, *.cpp)
	cd $(DOC_DIR) && $(DOXYGEN)

doc: $(DOC_DIR)/doxygen/html/index.html

clean:;		rm -Rf $(OBJ) $(LIBARCH) $(BINS) $(PYTHON_CLEAN) $(TEST_DIR)/test

