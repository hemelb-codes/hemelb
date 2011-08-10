include $(MK)/header.mk

TARGETS := CassonRheologyModel.o \
           TruncatedPowerLawRheologyModel.o
           
INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
