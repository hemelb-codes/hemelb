include $(MK)/header.mk

TARGETS := Implementation.o \
	   HFunction.o \
	   CollisionOperator.o \
	   LBGK.o \

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
