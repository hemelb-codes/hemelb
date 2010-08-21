# -*- mode: makefile; -*-
CXX := mpic++

EXE := hemelb

HEMELB_DEFS := CCS BSD DARWIN

HEMELB_CXXFLAGS := -O3 -pthread -Wunused

HEMELB_INCLUDEPATHS := /usr/local/include

HEMELB_LIBPATHS := /usr/local/lib 

# vim: set ft=make :