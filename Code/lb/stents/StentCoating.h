
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_STENTS_STENTCOATING_H
#define HEMELB_LB_STENTS_STENTCOATING_H

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
      class StentCoating : public Stent
      {
        public:
          StentCoating();
          virtual ~StentCoating();
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

          void SetDensity(const LatticeDensity& rho)
          {
            densityMean = rho;
          }

          distribn_t GetCoatingDiffusivity() const
          {
            return Dc;
          }

          distribn_t GetCoatingThickness() const
          {
            return lc;
          } 

          void SetCoatingDiffusivity(const distribn_t& D)
          {
            Dc = D;
          }

          void SetCoatingThickness(const distribn_t& l)
          {
            lc = l;
          }

        private:

          LatticeDensity densityMean;
          distribn_t Dc;
          distribn_t lc;
      };

    }
  }
}

#endif /* HEMELB_LB_STENTS_STENTCOATING_H */
