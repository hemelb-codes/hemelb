
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "lb/stents/StentConstant.h"
#include "configuration/SimConfig.h"
#include "net/IOCommunicator.h"

namespace hemelb
{
  namespace lb
  {
    namespace stents
    {

      StentConstant::StentConstant() :
        Stent(), densityMean(1.0)
      {

      }

      Stent* StentConstant::Clone() const
      {
        StentConstant* copy = new StentConstant(*this);

        return copy;
      }

      StentConstant::~StentConstant()
      {

      }

      LatticeDensity StentConstant::GetDensity(LatticeTimeStep time_step) const
      {
        // Calculate the target
        LatticeDensity target = densityMean;

        return target;
      }

    }
  }
}
