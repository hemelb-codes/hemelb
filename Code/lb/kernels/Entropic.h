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
      template<class LatticeType> class Entropic;

      template<class LatticeType>
      struct HydroVars<Entropic<LatticeType> > : public HydroVarsBase<LatticeType>
      {
        public:
          HydroVars(const distribn_t* const f) :
              HydroVarsBase<LatticeType>(f)
          {

          }

          site_t index;
      };

      /**
       * Entropic: This class implements the entropic kernel, a modification to the standard
       * LBGK kernel which ensures the increase of entropy.
       */
      template<class LatticeType>
      class Entropic : public BaseKernel<Entropic<LatticeType>, LatticeType>
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

          inline void DoCalculateDensityVelocityFeq(HydroVars<Entropic<LatticeType> >& hydroVars, site_t index)
          {
            hydroVars.index = index;
            LatticeType::CalculateEntropicDensityVelocityFEq(hydroVars.f,
                                                             hydroVars.density,
                                                             hydroVars.v_x,
                                                             hydroVars.v_y,
                                                             hydroVars.v_z,
                                                             hydroVars.f_eq.f);

            for (unsigned int ii = 0; ii < LatticeType::NUMVECTORS; ++ii)
            {
              hydroVars.f_neq.f[ii] = hydroVars.f[ii] - hydroVars.f_eq.f[ii];
            }
          }

          inline void DoCalculateFeq(HydroVars<Entropic<LatticeType> >& hydroVars, site_t index)
          {
            hydroVars.index = index;
            LatticeType::CalculateEntropicFeq(hydroVars.density,
                                              hydroVars.v_x,
                                              hydroVars.v_y,
                                              hydroVars.v_z,
                                              hydroVars.f_eq.f);

            for (unsigned int ii = 0; ii < LatticeType::NUMVECTORS; ++ii)
            {
              hydroVars.f_neq.f[ii] = hydroVars.f[ii] - hydroVars.f_eq.f[ii];
            }
          }

          inline void DoCollide(const LbmParameters* const lbmParams, HydroVars<Entropic<LatticeType> >& hydroVars)
          {
            distribn_t alpha = CalculateAlpha(lbmParams->GetTau(), hydroVars, oldAlpha[hydroVars.index]);
            oldAlpha[hydroVars.index] = alpha;

            for (Direction direction = 0; direction < LatticeType::NUMVECTORS; ++direction)
            {
              hydroVars.GetFPostCollision()[direction] = hydroVars.f[direction]
                  + (alpha * lbmParams->GetBeta()) * hydroVars.f_neq.f[direction];
            }
          }

          inline void Reset(InitParams* init)
          {
            for (site_t i = 0; i < init->siteCount; i++)
            {
              oldAlpha[i] = 2.0;
            }
          }

        private:
          distribn_t* oldAlpha;

          double CalculateAlpha(const distribn_t tau,
                                const HydroVars<Entropic<LatticeType> >& hydroVars,
                                double prevAlpha)
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

      };

    }
  }
}

#endif /* HEMELB_LB_KERNELS_ENTROPIC_H */
