#include "extraction/WholeGeometrySelector.h"

namespace hemelb
{
  namespace extraction
  {
    bool WholeGeometrySelector::IsWithinGeometry(const extraction::IterableDataSource& data,
                                                 const util::Vector3D<float>& location)
    {
      return true;
    }
  }
}
