include $(MK)/header.mk

TARGETS := Implementation.o \
	   CollisionOperator.o \
	   LBGK.o \
	   ELBM.o \

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
