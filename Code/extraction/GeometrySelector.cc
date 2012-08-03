// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

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
