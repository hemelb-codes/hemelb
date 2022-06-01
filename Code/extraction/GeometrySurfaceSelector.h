// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_EXTRACTION_GEOMETRYSURFACESELECTOR_H
#define HEMELB_EXTRACTION_GEOMETRYSURFACESELECTOR_H

#include "extraction/GeometrySelector.h"

namespace hemelb::extraction
{

  class GeometrySurfaceSelector : public GeometrySelector
  {
  public:
    ~GeometrySurfaceSelector() override = default;
    GeometrySelector* clone() const override;

  protected:
    bool IsWithinGeometry(const extraction::IterableDataSource& data, const util::Vector3D<site_t>& location) const override;
  };

}

#endif // HEMELB_EXTRACTION_GEOMETRYSURFACESELECTOR_H
