
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_KERNELS_ENTROPIC_H
#define HEMELB_LB_KERNELS_ENTROPIC_H

#include "lb/kernels/BaseKernel.h"
#include "lb/HFunction.h"
#include "util/utilityFunctions.h"

namespace hemelb
{
  namespace lb
  {
    namespace kernels
    {
      /**
       * Entropic: an incomplete kernel implementation that includes some common elements for implementations of
       * entropic LB.
       */
      template<typename LatticeType>
      class Entropic
      {
        public:
          /**
           * Destructor, deallocates the alpha array.
           */
          ~Entropic()
          {
            delete[] oldAlpha;
          }

          /**
           * Performs the Entropic LB collision (using alpha as a relaxation parameter)
           * @param lbmParams
           * @param hydroVars
           */
          template<typename HydroVarsType>
          inline void DoCollide(const LbmParameters* const lbmParams, HydroVarsType& hydroVars)
          {
            distribn_t alpha = CalculateAlpha(lbmParams->GetTau(), hydroVars, oldAlpha[hydroVars.index]);
            oldAlpha[hydroVars.index] = alpha;

            for (Direction direction = 0; direction < LatticeType::NUMVECTORS; ++direction)
            {
              hydroVars.SetFPostCollision(direction,
                                          hydroVars.f[direction]
                                              + (alpha * lbmParams->GetBeta()) * hydroVars.f_neq.f[direction]);
            }
          }

        protected:
          /**
           * Constructs the alpha array.
           * @param initParams
           */
          Entropic(InitParams* initParams)
          {
            oldAlpha = new distribn_t[initParams->latDat->GetLocalFluidSiteCount()];
            // Initialises the value of alpha to 2.0 for every site.
            for (site_t i = 0; i < initParams->latDat->GetLocalFluidSiteCount(); i++)
            {
              oldAlpha[i] = 2.0;
            }
          }

          /**
           * Calculates the new value of alpha (the relaxation parameter) and updates the cache.
           * @param tau The value of tau to use as a basis for calculating alpha.
           * @param hydroVars
           * @param prevAlpha Alpha value from the previous timestep.
           * @return
           */
          template<typename HydroVarsType>
          double CalculateAlpha(const distribn_t tau, const HydroVarsType& hydroVars, double prevAlpha)
          {
            bool big = false;
            double deviation = 0.0;

            for (unsigned int i = 0; i < LatticeType::NUMVECTORS; i++)
            {
              // Papers suggest f_eq - f < 0.001 or (f_eq - f)/f < 0.01 for the point to have approx alpha = 2
              // Accuracy can change depending on stability requirements, because the more NR evaluations it skips
              // the more of the simulation is in the LBGK limit.
              deviation = util::NumericalFunctions::max(fabs( (hydroVars.f_eq.f[i] - hydroVars.f[i]) / hydroVars.f[i]),
                                                        deviation);
              if (deviation > 1.0E-2)
              {
                big = true;
                // if this is the case Brent isn't called so the value of deviation doesn't matter
                break;
              }
            }

            if (big)
            {
              HFunction<LatticeType> HFunc(hydroVars.f, hydroVars.f_eq.f);

              // This is in case previous Alpha was calculated to be zero (does happen occasionally if f_eq - f is small
              prevAlpha = (prevAlpha < 2.0 * tau ?
                2.0 :
                prevAlpha);

              return (hemelb::util::NumericalMethods::NewtonRaphson(&HFunc, prevAlpha, 1.0E-6));
            }
            else
            {
              // Happens a lot near equilibrium. In this limit we return LBGK to avoid unnecessary calculations
              // Should really be a small number to accommodate round off and truncations errors, but only 0.0
              // causes problems as it causes a NaN when bracketing in Brent.
              if (deviation == 0.0)
              {
                return 2.0;
              }

              HFunction<LatticeType> HFunc(hydroVars.f, hydroVars.f_eq.f);

              // The bracket is very large, but it should guarantee that a root is enclosed
              double alphaLower = 2.0 * (tau), HLower;
              double alphaHigher = 2.0 * (tau) / deviation, HHigher;

              HFunc(alphaLower, HLower);
              HFunc(alphaHigher, HHigher);

              // The root should be enclosed, but in case it isn't return some default
              // Chosen to return 2.0 as that is the LBGK case
              // Very often if f is v close to equilibrium a root will not be enclosed (rounding and truncation errors)
              // Doesn't really matter what is returned then as f_neq is negligible in that case
              if (HLower * HHigher >= 0.0)
              {
                return 2.0;
              }

              return (hemelb::util::NumericalMethods::Brent(&HFunc,
                                                            alphaLower,
                                                            HLower,
                                                            alphaHigher,
                                                            HHigher,
                                                            1.0E-6,
                                                            1.0E-12));
            }

          }

          /**
           * Stores the value of alpha (the relaxation parameter) from the previous iteration.
           */
          distribn_t* oldAlpha;
      };
    }
  }
}

#endif /* HEMELB_LB_KERNELS_ENTROPIC_H */
