
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_STENTS_STENTCONSTANT_H
#define HEMELB_LB_STENTS_STENTCONSTANT_H

#include "lb/stents/Stent.h"

namespace hemelb
{
  namespace lb
  {
    namespace stents
    {

      /*
       * Template values are chosen to be tUpdatePeriod = 1 and tComms = false
       * The cosine pressure trace can be easily calculated locally by any proc and
       * there is no need to store values for any time step beyond the current one
       */
      class StentConstant : public Stent
      {
        public:
          StentConstant();
          virtual ~StentConstant();
          virtual Stent* Clone() const;
          virtual void Reset(SimulationState &state)
          {
            //pass;
          }

          LatticeDensity GetDensity(unsigned long time_step) const;

          LatticeDensity GetDensityMin() const
          {
            return densityMean;
          }
          LatticeDensity GetDensityMax() const
          {
            return densityMean;
          }

        private:

          LatticeDensity densityMean;
      };

    }
  }
}

#endif /* HEMELB_LB_STENTS_STENTCONSTANT_H */
