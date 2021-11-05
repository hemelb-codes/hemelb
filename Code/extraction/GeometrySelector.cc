// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "extraction/GeometrySelector.h"

namespace hemelb
{
  namespace extraction
  {
    bool GeometrySelector::Include(const extraction::IterableDataSource& data,
                                   const util::Vector3D<site_t>& location) const
    {
      if (!data.IsValidLatticeSite(location) || !data.IsAvailable(location))
      {
        return false;
      }

      return IsWithinGeometry(data, location);
    }

    util::Vector3D<float> GeometrySelector::LatticeToPhysical(const extraction::IterableDataSource& data, const util::Vector3D<site_t>& location) const {
      return util::Vector3D<float>{location} * float(data.GetVoxelSize()) + data.GetOrigin().as<float>();
    }
  }
}
