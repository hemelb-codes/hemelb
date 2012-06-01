#include "GeometrySurfaceSelector.h"
#include <iostream>

namespace hemelb
{
  namespace extraction
  {
    bool GeometrySurfaceSelector::IsWithinGeometry(const extraction::IterableDataSource& data,
                                                   const util::Vector3D<site_t>& location)
    {
      return data.IsBoundarySite(location);
    }

  } /* namespace extraction */
} /* namespace hemelb */
