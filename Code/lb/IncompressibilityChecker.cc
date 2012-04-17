#include "lb/lattices/D3Q15.h"
#include "lb/IncompressibilityChecker.hpp"

namespace hemelb
{
  namespace lb
  {
    // Explicit instantiation
    template class IncompressibilityChecker<net::PhasedBroadcastRegular<> > ;
  }
}
