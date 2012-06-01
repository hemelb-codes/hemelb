#include "extraction/GeometrySelector.h"

namespace hemelb
{
  namespace extraction
  {
    bool GeometrySelector::Include(const extraction::IterableDataSource& data,
                                   const util::Vector3D<site_t>& location)
    {
      if (!data.IsValidLatticeSite(location) || !data.IsAvailable(location))
      {
        return false;
      }

      return IsWithinGeometry(data, location);
    }
  }
}
