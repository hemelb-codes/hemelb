
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "lb/stents/StentCoating.h"
#include "configuration/SimConfig.h"
#include "net/IOCommunicator.h"

namespace hemelb
{
  namespace lb
  {
    namespace stents
    {

      StentCoating::StentCoating() :
        Stent(), densityMean(1.0), Dc(1.0), lc(1.0)
      {

      }

      Stent* StentCoating::Clone() const
      {
        StentCoating* copy = new StentCoating(*this);

        return copy;
      }

      StentCoating::~StentCoating()
      {

      }

      LatticeDensity StentCoating::GetDensity(LatticeTimeStep time_step) const
      {
        // Calculate the target
        LatticeDensity target = densityMean;

        return target;
      }

    }
  }
}
