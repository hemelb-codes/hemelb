
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#include "lb/stents/StentConstantFlux.h"
#include "configuration/SimConfig.h"

namespace hemelb
{
  namespace lb
  {
    namespace stents
    {
      StentConstantFlux::StentConstantFlux() :
          maxSpeed(0.)
      {
      }

      StentConstantFlux::~StentConstantFlux()
      {
      }

      Stent* StentConstantFlux::Clone() const
      {
        Stent* copy = new StentConstantFlux(*this);
        return copy;
      }

      LatticeSpeed StentConstantFlux::GetFlux(const LatticeTimeStep t) const
      {
        LatticeSpeed max = maxSpeed;

        return max;
      }
    }
  }
}
