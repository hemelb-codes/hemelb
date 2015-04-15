// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_EXTRACTION_WHOLEGEOMETRYSELECTOR_H
#define HEMELB_EXTRACTION_WHOLEGEOMETRYSELECTOR_H

#include "extraction/GeometrySelector.h"

namespace hemelb
{
  namespace extraction
  {
    /**
     * Class for selecting the whole geometry.
     */
    class WholeGeometrySelector : public GeometrySelector
    {
      protected:
        /**
         * Returns true for all locations.
         *
         * @param data
         * @param location
         * @return
         */
        bool IsWithinGeometry(const extraction::IterableDataSource& data,
                              const util::Vector3D<site_t>& location);
    };
  }
}

#endif /* HEMELB_EXTRACTION_WHOLEGEOMETRYSELECTOR_H */
