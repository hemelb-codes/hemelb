#include "reporting/Timers.hpp"
namespace hemelb
{
  namespace reporting
  {
    template class TimersBase<HemeLBClockPolicy, MPICommsPolicy>; // explicit instantiate
  }
}
