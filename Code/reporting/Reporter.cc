#include "Reporter.hpp"
namespace hemelb
{
  namespace reporting
  {
    template class ReporterBase<Timers,FileWriterPolicy,MPICommsPolicy>; // explicit instantiate
  }
}
