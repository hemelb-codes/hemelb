# -*- mode: makefile; -*-
include $(MK)/config-default.mk

CXX := mpic++.openmpi

EXE := hemelb

HEMELB_DEFS := CCS

HEMELB_CXXFLAGS :=-O4 -pthread -Wunused -g

# vim: set ft=make :
