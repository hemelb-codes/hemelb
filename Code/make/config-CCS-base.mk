include $(MK)/config-linux.mk

CXX := mpic++.openmpi

HEMELB_DEFS += CCS

HEMELB_CXXFLAGS += -pthread -Wno-long-long -Wno-unused-result