include $(MK)/config-linux.mk

CXX := CC
HEMELB_CXXFLAGS := -g -fastsse --display_error_number --diag_suppress 611

HEMELB_STEERING_LIB := none
HEMELB_DEFS += NO_STREAKLINES

#$(EXE)_LIBS = parmetis/ParMetis-3.1.1/libparmetis.a parmetis/ParMetis-3.1.1/libmetis.a
#HEMELB_INCLUDEPATHS += parmetis/ParMetis-3.1.1
