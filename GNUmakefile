# Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
# SPDX-License-Identifier: BSD-2-Clause

# Portable GNU make build for Linux, macOS, and BSD.

SHELL := /bin/sh

CXX ?= c++
AR ?= ar
ARFLAGS := rcs
DOXYGEN ?= doxygen
UV ?= uv
CATCH_DIR := contrib/Catch2
CATCH_CPP := $(CATCH_DIR)/extras/catch_amalgamated.cpp

CPPFLAGS += -Icontrib/eigen -Isrc/libpdb++ -Isrc/zdock -Isrc/common \
            -Isrc/pdb -Iinclude
CXXFLAGS ?= -O3
CXXSTANDARD ?= -std=c++14
WARNINGS ?= -Wall -Wextra -pedantic
LDFLAGS ?=
LDLIBS ?=

BIN_DIR := bin
OBJ_DIR := build
LIB_DIR := lib
TEST_DIR := test
PYTHON_DIR := python

LIBRARY := zdock
LIBARCH := $(LIB_DIR)/lib$(LIBRARY).a

LIB_SOURCES := \
  src/libpdb++/pdbinput.cpp \
  src/libpdb++/pdb_read.cpp \
  src/libpdb++/pdb++.cpp \
  src/libpdb++/pdb_sscanf.cpp \
  src/libpdb++/pdb_type.cpp \
  src/libpdb++/pdb_sprntf.cpp \
  src/libpdb++/pdb_chars.cpp \
  src/zdock/TransformMultimer.cpp \
  src/zdock/Constraints.cpp \
  src/zdock/TransformLigand.cpp \
  src/zdock/TransformUtil.cpp \
  src/zdock/ZDOCK.cpp \
  src/pdb/PDB.cpp

TOOL_NAMES := createlig createmultimer pruning constraints centroids zdsplit zdunsplit
TOOL_SOURCES := \
  src/CreateLigand.cpp \
  src/CreateMultimer.cpp \
  src/Pruning.cpp \
  src/FilterConstraints.cpp \
  src/Centroids.cpp \
  src/Split.cpp \
  src/UnSplit.cpp
TEST_SOURCES := test/Test.cpp test/constraints.cpp test/formatting.cpp \
  test/multimer.cpp test/pdb.cpp test/prediction.cpp test/rotation.cpp \
  test/utils.cpp test/zdock.cpp test/zdock_errors.cpp

LIB_OBJECTS := $(patsubst %.cpp,$(OBJ_DIR)/%.o,$(LIB_SOURCES))
TOOL_OBJECTS := $(patsubst %.cpp,$(OBJ_DIR)/%.o,$(TOOL_SOURCES))
TEST_OBJECTS := $(patsubst %.cpp,$(OBJ_DIR)/%.o,$(TEST_SOURCES))
OBJECTS := $(LIB_OBJECTS) $(TOOL_OBJECTS) $(TEST_OBJECTS)
DEPS := $(OBJECTS:.o=.d)
BINS := $(addprefix $(BIN_DIR)/,$(TOOL_NAMES))
TEST_BIN := $(TEST_DIR)/test

.PHONY: all check-deps test cpp-test python-check python-test doc clean

all: check-deps $(BINS)

check-deps:
	@test -f contrib/eigen/Eigen/Dense || { \
	  echo "Eigen submodule missing; run: git submodule update --init --recursive" >&2; \
	  exit 1; \
	}
	@test -f $(CATCH_CPP) || { \
	  echo "Catch2 submodule missing; run: git submodule update --init --recursive" >&2; \
	  exit 1; \
	}

$(LIBARCH): $(LIB_OBJECTS) | $(LIB_DIR)
	$(AR) $(ARFLAGS) $@ $^

$(OBJ_DIR)/test/%.o: CPPFLAGS += -I$(CATCH_DIR)/extras \
  -DDATADIR=$(abspath $(TEST_DIR)/data)

$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(CXXSTANDARD) $(WARNINGS) -MMD -MP -c -o $@ $<

define LINK_TOOL
$(BIN_DIR)/$(1): $(OBJ_DIR)/$(2:.cpp=.o) $(LIBARCH) | $(BIN_DIR)
	$$(CXX) $$(LDFLAGS) -o $$@ $$^ $$(LDLIBS)
endef

$(eval $(call LINK_TOOL,createlig,src/CreateLigand.cpp))
$(eval $(call LINK_TOOL,createmultimer,src/CreateMultimer.cpp))
$(eval $(call LINK_TOOL,pruning,src/Pruning.cpp))
$(eval $(call LINK_TOOL,constraints,src/FilterConstraints.cpp))
$(eval $(call LINK_TOOL,centroids,src/Centroids.cpp))
$(eval $(call LINK_TOOL,zdsplit,src/Split.cpp))
$(eval $(call LINK_TOOL,zdunsplit,src/UnSplit.cpp))

$(TEST_BIN): $(TEST_OBJECTS) $(LIBARCH) $(OBJ_DIR)/catch_amalgamated.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

$(OBJ_DIR)/catch_amalgamated.o: $(CATCH_CPP)
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(CXXSTANDARD) $(WARNINGS) \
	  -I$(CATCH_DIR)/extras -c -o $@ $<

cpp-test: check-deps $(TEST_BIN)
	$(TEST_BIN)

python-test: all
	$(UV) --directory $(PYTHON_DIR) run python -m unittest discover -s tests -v

python-check:
	$(UV) --directory $(PYTHON_DIR) run isort --check-only src tests
	$(UV) --directory $(PYTHON_DIR) run black --check src tests
	$(UV) --directory $(PYTHON_DIR) run mypy src tests
	PYLINTHOME=$(abspath $(PYTHON_DIR)/.pylint.d) \
		$(UV) --directory $(PYTHON_DIR) run pylint src tests

test: cpp-test python-test

doc:
	cd doc && $(DOXYGEN)

$(BIN_DIR) $(LIB_DIR):
	mkdir -p $@

clean:
	rm -rf $(OBJ_DIR) $(LIB_DIR) $(BIN_DIR) $(TEST_BIN) \
	  $(PYTHON_DIR)/build $(PYTHON_DIR)/dist $(PYTHON_DIR)/src/*.egg-info \
	  $(PYTHON_DIR)/.mypy_cache $(PYTHON_DIR)/.pylint.d \
	  $(PYTHON_DIR)/src/zdock/__pycache__ $(PYTHON_DIR)/tests/__pycache__

-include $(DEPS)
