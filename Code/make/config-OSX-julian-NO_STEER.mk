HEMELB_STEERING_LIB := none

include $(MK)/config-OSX-base.mk
HEMELB_CXXFLAGS := -g -pedantic -Wall -Wextra -Wno-unused -O2

include $(MK)/config-OSX-julian-base.mk
