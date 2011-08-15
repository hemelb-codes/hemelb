include $(MK)/header.mk

TARGETS := LBGK.o \
	   ELBM.o \

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
