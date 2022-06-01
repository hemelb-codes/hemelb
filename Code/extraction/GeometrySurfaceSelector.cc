// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "extraction/GeometrySurfaceSelector.h"

namespace hemelb::extraction
{

  GeometrySelector* GeometrySurfaceSelector::clone() const {
    return new GeometrySurfaceSelector{};
  }

  bool GeometrySurfaceSelector::IsWithinGeometry(const extraction::IterableDataSource& data,
						 const util::Vector3D<site_t>& location) const
  {
    return data.IsWallSite(location);
  }

}
