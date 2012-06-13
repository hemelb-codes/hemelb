add_library(hemelb_net IteratedAction.cc BaseNet.cc 
mixins/CoalescePointPoint.cc
mixins/SeparatedGathers.cc 
mixins/GathersViaPointPoint.cc
mixins/SeparatedAllToAll.cc
mixins/StoringNet.cc ProcComms.cc
phased/StepManager.cc)
