include $(MK)/config-linux.mk

CXX := mpic++.openmpi
CC := mpicc.openmpi

HEMELB_DEFS += CCS

HEMELB_CXXFLAGS += -pthread -Wno-long-long -Wno-unused-result -Wconversion -Wall -Wextra
HEMELB_CCFLAGS += -pthread -Wno-long-long -Wno-unused-result -Wconversion

HEMELB_LOG_LEVEL = warning
