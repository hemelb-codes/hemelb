include $(MK)/config-linux.mk
# Note that as of Nov 2010, you need to use the intel and mvapich2 modules to compile
CXX := mpicxx
AR := xiar

HEMELB_CXXFLAGS := -O3 -xW -ipo -ipo-jobs4 -pthread

# The MPICH define is to get around a bug in the MPI2 C++ binding on MPICH
HEMELB_DEFS += MPICH_IGNORE_CXX_SEEK


