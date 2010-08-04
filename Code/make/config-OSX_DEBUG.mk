# -*- mode: makefile; -*-
CXX := /usr/local/bin/mpic++

EXE := hemelb

HEMELB_DEFS := CCS BSD DARWIN DEBUG

HEMELB_CXXFLAGS := -g -pthread -Wunused

HEMELB_INCLUDEPATHS := /usr/local/include

HEMELB_LIBPATHS := /usr/local/lib 

# vim: set ft=make :