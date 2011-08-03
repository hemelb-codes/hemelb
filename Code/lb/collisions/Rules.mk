include $(MK)/header.mk

SUBDIRS := implementations

TARGETS := CollisionVisitor.o \
           CollisionOperator.o

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
