SHELL		= /bin/sh
PREFIX  = $(HOME)/local
CXX = clang++
STRIP = strip
BIN = bin/test
SRC = -I/opt/local/include/eigen3
BINOBJ = src/test.cpp
OBJ_DIR := build
OPT		= -O3 -DEIGEN_USE_LAPACKE
DEBUG		=
CXXFLAGS		= $(OPT) $(DEBUG)
LIBRARY		= pdb++
LIBARCH		= lib/lib$(LIBRARY).a
BIN_SOURCES = src/test.cpp
LIB_SOURCES = src/PDB.cpp src/pdb++.cpp src/pdb_chars.cpp src/pdb_read.cpp src/pdb_sprntf.cpp src/pdb_sscanf.cpp src/pdb_type.cpp src/pdbinput.cpp

rwildcard=$(foreach d,$(wildcard $1*),$(call rwildcard,$d/,$2) $(filter $(subst *,%,$2),$d))
HEADERS = $(call rwildcard, src/, *.h) $(call rwildcard, src/, *.hpp) $(call rwildcard, src/, *.i)
OBJ := $(addprefix $(OBJ_DIR)/, $(patsubst %.cpp, %.o, $(call rwildcard, src/, *.cpp)))
BINOBJ = $(addprefix $(OBJ_DIR)/, $(patsubst %.cpp, %.o, $(BIN_SOURCES)))
LIBOBJ = $(addprefix $(OBJ_DIR)/, $(patsubst %.cpp, %.o, $(LIB_SOURCES)))

all:		$(OBJ_DIR) $(BIN) Makefile

$(LIBARCH): $(LIBOBJ)
	ar cru $(LIBARCH) $(LIBOBJ)

$(OBJ_DIR)/%.o: %.cpp $(HEADERS)
	@mkdir -p $(OBJ_DIR)/$(shell dirname $<)
	$(CXX) $(CXXFLAGS) $(SRC) -c -o $@ $<

$(BIN): $(BINOBJ) $(LIBARCH)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LD_FLAGS)
	$(STRIP) $@

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

clean:;		rm -f $(OBJ) $(LIBARCH) $(BIN)

