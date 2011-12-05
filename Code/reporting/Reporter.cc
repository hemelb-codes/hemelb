#include "Reporter.hpp"
namespace hemelb
{
  namespace reporting
  {
    template class ReporterBase<Timers, FileWriterPolicy, MPICommsPolicy, net::PhasedBroadcastRegular<> > ; // explicit instantiate
  }
}
