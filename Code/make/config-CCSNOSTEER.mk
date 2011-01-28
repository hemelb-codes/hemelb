# -*- mode: makefile; -*-
CXX := mpic++.openmpi

EXE := hemelb

HEMELB_DEFS := CCS NO_STEER

HEMELB_CXXFLAGS :=-O3 -pthread -Wunused

# vim: set ft=make :
