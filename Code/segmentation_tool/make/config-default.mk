CXX := g++

EXE := segtool

($EXE)_DEBUG_LEVEL := 0
HEMELB_CXXFLAGS := -g -pedantic -Wall -Wextra -Wno-unused
HEMELB_DEFS :=

$(EXE)_LIBS := -lGLU -lGL -lglut
