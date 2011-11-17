include $(MK)/config-linux.mk

CXX := mpiCC
EXE := hemelb

HEMELB_CXXFLAGS := -O4
HEMELB_DEFS += BENCH NOOPENMP
