#include "Decomposition.hpp"

namespace hemelb
{
  namespace geometry
  {
    template class DecompositionBase<net::Net, topology::NetworkTopology>; // explicit instantiate
  }
}