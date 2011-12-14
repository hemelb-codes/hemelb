#include "lb/lattices/Lattice.h"

namespace hemelb
{
  namespace lb
  {
    namespace lattices
    {
      template<typename T>
      LatticeInfo* Lattice<T>::singletonInfo = NULL;
    }
  }
}
