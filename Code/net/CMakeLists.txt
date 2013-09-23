add_library(hemelb_net mpi.cc IteratedAction.cc BaseNet.cc 
mixins/pointpoint/CoalescePointPoint.cc
mixins/pointpoint/SeparatedPointPoint.cc
mixins/pointpoint/ImmediatePointPoint.cc
mixins/gathers/SeparatedGathers.cc 
mixins/gathers/ViaPointPointGathers.cc
mixins/alltoall/SeparatedAllToAll.cc
mixins/alltoall/ViaPointPointAllToAll.cc
mixins/StoringNet.cc ProcComms.cc
phased/StepManager.cc)
configure_file (
  "${PROJECT_SOURCE_DIR}/net/BuildInfo.h.in"
  "${PROJECT_BINARY_DIR}/net/BuildInfo.h"
  )