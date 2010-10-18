# -*- mode: makefile; -*-
include $(MK)/config-default.mk

CXX := mpic++.openmpi

EXE := hemelb

HEMELB_DEBUG_LEVEL := 1
HEMELB_DEFS := CCS

HEMELB_CXXFLAGS :=-O4 -pthread -Wunused -g

# vim: set ft=make :
