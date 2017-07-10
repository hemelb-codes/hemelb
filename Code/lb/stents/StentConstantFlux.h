
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_STENTS_STENTCONSTANTFLUX_H
#define HEMELB_LB_STENTS_STENTCONSTANTFLUX_H
#include "lb/stents/StentFlux.h"

namespace hemelb
{
  namespace lb
  {
    namespace stents
    {
      class StentConstantFlux : public StentFlux
      {
        public:
          StentConstantFlux();
          virtual ~StentConstantFlux();
          Stent* Clone() const;
          LatticeSpeed GetFlux(const LatticeTimeStep t) const;

          void SetFlux(const LatticeSpeed& flux)
          {
            maxSpeed = flux;
          }

        protected:
          LatticeSpeed maxSpeed;
      };
    }
  }
}
#endif // HEMELB_LB_STENTS_STENTCONSTANTFLUX_H
