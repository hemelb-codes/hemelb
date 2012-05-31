add_library(hemelb_net IteratedAction.cc net.cc BaseNet.cc 
mixins/CoalescePointPoint.cc
mixins/SeparatedGathers.cc 
mixins/GathersViaPointPoint.cc
mixins/StoringNet.cc ProcComms.cc)
