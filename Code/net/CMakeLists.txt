add_library(hemelb_net
  MpiDataType.cc MpiEnvironment.cc MpiError.cc
  MpiCommunicator.cc MpiGroup.cc MpiFile.cc
 IteratedAction.cc BaseNet.cc 
IOCommunicator.cc
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
configure_file (
  "${PROJECT_SOURCE_DIR}/net/MpiConstness.h.in"
  "${PROJECT_BINARY_DIR}/net/MpiConstness.h"
  )