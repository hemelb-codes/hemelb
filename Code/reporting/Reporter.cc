#include "Reporter.hpp"
namespace hemelb
{
  namespace reporting
  {
    template class ReporterBase<HemeLBClockPolicy, FileWriterPolicy, MPICommsPolicy, net::PhasedBroadcastRegular<> > ; // explicit instantiate
  }
}
