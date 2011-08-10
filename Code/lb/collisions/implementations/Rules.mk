include $(MK)/header.mk

TARGETS := LBGK.o \
           ELBM.o \
           LBGKNN.o

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
