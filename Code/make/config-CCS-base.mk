include $(MK)/config-linux.mk

CXX := mpic++.openmpi
CC := mpicc.openmpi

HEMELB_DEFS += CCS

HEMELB_CXXFLAGS += -pthread -Wno-long-long -Wno-unused-result
