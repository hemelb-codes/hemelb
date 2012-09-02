add_library(hemelb_net IteratedAction.cc BaseNet.cc 
mixins/CoalescePointPoint.cc
mixins/SeparatedPointPoint.cc
mixins/SeparatedGathers.cc 
mixins/GathersViaPointPoint.cc
mixins/SeparatedAllToAll.cc
mixins/AllToAllViaPointPoint.cc
mixins/StoringNet.cc ProcComms.cc
phased/StepManager.cc)
