#ifndef HEMELB_LB_KERNELS_ENTROPIC_H
#define HEMELB_LB_KERNELS_ENTROPIC_H

#include <cstdlib>

#include "lb/kernels/BaseKernel.h"
#include "lb/HFunction.h"
#include "util/utilityFunctions.h"

namespace hemelb
{
  namespace lb
  {
    namespace kernels
    {
      // We have to declare this up here in order for it to be used as a template parameter in the
      // following declaration. Moving the template specialisation to the bottom of the file would
      // prevent it from being used as the HydroVars for this kernel.
      class Entropic;

      template<>
      struct HydroVars<Entropic> : public HydroVarsBase
      {
        public:
          HydroVars(const distribn_t* const f) :
            HydroVarsBase(f)
          {

          }

          distribn_t alpha;
          site_t index;
      };

      /**
       * Entropic: This class implements the entropic kernel, a modification to the standard
       * LBGK kernel which ensures the increase of entropy.
       */
      class Entropic : public BaseKernel<Entropic>
      {
        public:
          Entropic(InitParams& initParams)
          {
            oldAlpha = new distribn_t[initParams.siteCount];

            Reset(&initParams);
          }

          ~Entropic()
          {
            delete[] oldAlpha;
          }

          void DoCalculateDensityVelocityFeq(HydroVars<Entropic>& hydroVars, site_t index)
          {
            hydroVars.index = index;
            D3Q15::CalculateEntropicDensityVelocityFEq(hydroVars.f,
                                                       hydroVars.density,
                                                       hydroVars.v_x,
                                                       hydroVars.v_y,
                                                       hydroVars.v_z,
                                                       hydroVars.f_eq);
          }

          void DoCalculateFeq(HydroVars<Entropic>& hydroVars, site_t index)
          {
            hydroVars.index = index;
            D3Q15::CalculateEntropicFeq(hydroVars.density,
                                        hydroVars.v_x,
                                        hydroVars.v_y,
                                        hydroVars.v_z,
                                        hydroVars.f_eq);
          }

          distribn_t DoCollide(const LbmParameters* const lbmParams,
                               HydroVars<Entropic>& hydroVars,
                               unsigned int direction)
          {
            // TODO this calculation is being done in the wrong place and hence happens 15* too often.
            hydroVars.alpha
                = CalculateAlpha(lbmParams->Tau(), hydroVars, oldAlpha[hydroVars.index]);
            oldAlpha[hydroVars.index] = hydroVars.alpha;

            return hydroVars.f[direction] + (hydroVars.alpha * lbmParams->Beta())
                * hydroVars.f_neq[direction];
          }

          void Reset(InitParams* init)
          {
            for (site_t i = 0; i < init->siteCount; i++)
            {
              oldAlpha[i] = 2.0;
            }
          }

        private:
          distribn_t* oldAlpha;

          double CalculateAlpha(const distribn_t tau,
                                const HydroVars<Entropic>& hydroVars,
                                double prevAlpha)
          {
            bool big = false;
            double deviation = 0.0;

            for (unsigned int i = 0; i < D3Q15::NUMVECTORS; i++)
            {
              // Papers suggest f_eq - f < 0.001 or (f_eq - f)/f < 0.01 for the point to have approx alpha = 2
              // Accuracy can change depending on stability requirements, because the more NR evaluations it skips
              // the more of the simulation is in the LBGK limit.
              deviation = util::NumericalFunctions::max(fabs( (hydroVars.f_eq[i] - hydroVars.f[i])
                  / hydroVars.f[i]), deviation);
              if (deviation > 1.0E-2)
              {
                big = true;
                // if this is the case Brent isn't called so the value of deviation doesn't matter
                break;
              }
            }

            if (big)
            {
              HFunction HFunc(hydroVars.f, hydroVars.f_eq);

              // This is in case previous Alpha was calculated to be zero (does happen occasionally if f_eq - f is small
              prevAlpha = (prevAlpha < 2.0 * tau
                ? 2.0
                : prevAlpha);

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

              HFunction HFunc(hydroVars.f, hydroVars.f_eq);

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
                                                            1.0E-6));
            }

          }

      };

    }
  }
}

#endif /* HEMELB_LB_KERNELS_ENTROPIC_H */
