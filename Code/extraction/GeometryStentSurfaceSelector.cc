
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "extraction/GeometryStentSurfaceSelector.h"

namespace hemelb
{
  namespace extraction
  {
    bool GeometryStentSurfaceSelector::IsWithinGeometry(const extraction::IterableDataSource& data,
                                                   const util::Vector3D<site_t>& location)
    {
      return data.IsStentWallSite(location);
    }

  } /* namespace extraction */
} /* namespace hemelb */
