
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#include "lb/stents/StentFlux.h"
#include "configuration/SimConfig.h"

namespace hemelb
{
  namespace lb
  {
    namespace stents
    {

      StentFlux::StentFlux()
      {
      }

      StentFlux::~StentFlux()
      {
      }

      LatticeDensity StentFlux::GetDensityMin() const
      {
        return 1.0;
      }
      LatticeDensity StentFlux::GetDensityMax() const
      {
        return 1.0;
      }
      LatticeDensity StentFlux::GetDensity(LatticeTimeStep time_step) const
      {
        return 1.0;
      }
    }
  }
}
