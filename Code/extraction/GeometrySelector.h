// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_EXTRACTION_GEOMETRYSELECTOR_H
#define HEMELB_EXTRACTION_GEOMETRYSELECTOR_H

#include "units.h"
#include "extraction/IterableDataSource.h"
#include "util/Vector3D.h"

namespace hemelb
{
  namespace extraction
  {
    /**
     * Class representing a general geometry selector, i.e. a means of consistently selecting certain
     * lattice sites from a geometry.
     *
     * This is done with the Include function, overridden by implementors of the interface.
     */
    class GeometrySelector
    {
      public:
        /**
         * Virtual destructor.
         */
        virtual ~GeometrySelector()
        {

        }

        /**
         * Returns true if the given location is within the selection.
         */
        bool Include(const extraction::IterableDataSource& data, const util::Vector3D<site_t>& location) const;

      protected:
        virtual bool
            IsWithinGeometry(const extraction::IterableDataSource& data, const util::Vector3D<site_t>& location) const = 0;

        // Helper to do the coordinate transform for subclasses
        util::Vector3D<float> LatticeToPhysical(const extraction::IterableDataSource& data, const util::Vector3D<site_t>& location) const;
    };
  }
}

#endif /* HEMELB_EXTRACTION_GEOMETRYSELECTOR_H */
