include $(MK)/config-linux.mk

CXX := mpic++.openmpi
CC := mpicc.openmpi

HEMELB_DEFS += CCS

HEMELB_CXXFLAGS += -pthread -Wno-long-long -Wno-unused-result -Wconversion
HEMELB_CCFLAGS += -pthread -Wno-long-long -Wno-unused-result -Wconversion

HEMELB_LOG_LEVEL = warning

HEMELB_CXXFLAGS += -O4
HEMELB_CCFLAGS += -O4

PMETIS_LIBRARY_DIR := $(TOP)/parmetis/build/Linux-i686/libmetis/ $(TOP)/parmetis/build/Linux-i686/libparmetis/

