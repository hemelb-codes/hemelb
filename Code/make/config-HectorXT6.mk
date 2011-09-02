include $(MK)/config-linux.mk

CXX := CC
HEMELB_CXXFLAGS := -Mnounderscoring --one_instantiation_per_object --prelink_objects
#-g -fastsse --display_error_number --diag_suppress 611

HEMELB_STEERING_LIB := none
HEMELB_DEFS += NO_STREAKLINES

HEMELB_LOG_LEVEL := info