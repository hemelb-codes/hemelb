include $(MK)/config-linux.mk

CXX := mpiCC
EXE := hemelb
TARGETS := $(EXE)
UNITTESTS := 
NOTESTS := true
HEMELB_CXXFLAGS := -O4
HEMELB_DEFS += BENCH NOOPENMP
