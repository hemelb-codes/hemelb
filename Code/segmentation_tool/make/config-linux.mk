include $(MK)/config-default.mk

CXX := g++

EXE := segtool

HEMELB_DEBUG_LEVEL := 0
HEMELB_CXXFLAGS := -g -pedantic -Wall -Wextra -Wno-unused
HEMELB_DEFS := MESH

$(EXE)_DEBUG_LEVEL := 0
$(EXE)_CXXFLAGS := -g -pedantic -Wall -Wextra -Wno-unused
$(EXE)_DEFS := MESH

