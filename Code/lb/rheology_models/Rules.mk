include $(MK)/header.mk

TARGETS := AbstractRheologyModel.o \
		   CassonRheologyModel.o \
           TruncatedPowerLawRheologyModel.o \
           CarreauYasudaRheologyModel.o
           
INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
