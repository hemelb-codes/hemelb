#include "Reporter.hpp"
namespace hemelb
{
  namespace reporting
  {
    template class ReporterBase<HemeLBClockPolicy, MPICommsPolicy, net::PhasedBroadcastRegular<> > ; // explicit instantiate
  }
}
