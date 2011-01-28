include $(MK)/header.mk

SRCS := common.cc

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
