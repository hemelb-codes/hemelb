add_library(hemelb_net IteratedAction.cc BaseNet.cc 
mixins/CoalescePointPoint.cc
mixins/SeparatedGathers.cc 
mixins/GathersViaPointPoint.cc
mixins/SeparatedAllToAll.cc
mixins/AllToAllViaPointPoint.cc
mixins/StoringNet.cc ProcComms.cc
phased/StepManager.cc)
configure_file (
  "${PROJECT_SOURCE_DIR}/net/BuildInfo.h.in"
  "${PROJECT_BINARY_DIR}/net/BuildInfo.h"
  )