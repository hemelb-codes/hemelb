# Non-recursive makefile from
# http://code.google.com/p/nonrec-make/
# version 0.3

SHELL := /bin/bash
RUNDIR := $(shell pwd)
ifndef TOP
TOP := $(shell rd=$(RUNDIR); top=$$rd; \
               until [ -r $$top/Rules.top ]; do \
                 oldtop=$$top; \
                 cd ..; top=`pwd`; \
                 if [ $$oldtop = $$top ]; then \
                   top=$$rd; break; \
                 fi; \
               done; \
               echo $$top)
endif

MK := $(TOP)/make

.PHONY: all clean clean_all clean_tree

.DEFAULT_GOAL := dir_$(RUNDIR)
clean : clean_$(RUNDIR)
clean_tree : clean_tree_$(RUNDIR)

include $(TOP)/Rules.top
