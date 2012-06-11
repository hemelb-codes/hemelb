#include "extraction/GeometrySurfaceSelector.h"

namespace hemelb
{
  namespace extraction
  {
    bool GeometrySurfaceSelector::IsWithinGeometry(const extraction::IterableDataSource& data,
                                                   const util::Vector3D<site_t>& location)
    {
      return data.IsEdgeSite(location);
    }

  } /* namespace extraction */
} /* namespace hemelb */
