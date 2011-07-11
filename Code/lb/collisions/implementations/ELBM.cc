#include "lb/collisions/implementations/ELBM.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {
      namespace implementations
      {

        double* ELBM::alpha;
        size_t ELBM::currentAlphaIndex;

        void ELBM::createAlphaArray(const size_t size)
        {
          alpha = new double[size];
          for (size_t i = 0; i < size; i++)
            alpha[i] = 2.0;
        }

        void ELBM::getSiteValues(const distribn_t* f,
                                 distribn_t &density,
                                 distribn_t &v_x,
                                 distribn_t &v_y,
                                 distribn_t &v_z,
                                 distribn_t* f_eq,
                                 const site_t index)
        {
          currentAlphaIndex = index;
          D3Q15::CalculateEntropicDensityVelocityFEq(f, density, v_x, v_y, v_z, f_eq);
          alpha[index] = getAlpha(f, f_eq, alpha[index]);
        }

        void ELBM::getBoundarySiteValues(const distribn_t* f,
                                         const distribn_t &density,
                                         const distribn_t &v_x,
                                         const distribn_t &v_y,
                                         const distribn_t &v_z,
                                         distribn_t* f_eq,
                                         const site_t index)
        {
          currentAlphaIndex = index;
          D3Q15::CalculateEntropicFeq(density, v_x, v_y, v_z, f_eq);
          alpha[index] = getAlpha(f, f_eq, alpha[index]);
        }

        // Also updates lFEq_i to be lFNeq_i
        distribn_t ELBM::getOperatorElement(distribn_t &f_i,
                                            distribn_t &f_eq_i,
                                            const LbmParameters* iLbmParams)
        {
          return (alpha[currentAlphaIndex] * iLbmParams->Beta * (f_eq_i = f_i - f_eq_i));
        }

        double ELBM::getAlpha(const distribn_t* lF, const distribn_t* lFEq, double prevAlpha)
        {
          bool small = true;

          for (unsigned int i = 0; i < D3Q15::NUMVECTORS; i++)
          {
            // Papers suggest f_eq - f < 0.001 or (f_eq - f)/f < 0.01 for the point to have approx alpha = 2
            // Accuracy can change depending on stability requirements, because the more NR evaluations it skips
            // the more of the simulation is in the LBGK limit.
            if (fabs( (lFEq[i] - lF[i]) / lF[i]) > 1.0E-2)
            {
              small = false;
              break;
            }
          }

          HFunction HFunc(lF, lFEq);

          if (small)
          {
            return 2.0;

            double alphaLower = 1.8, HLower;
            double alphaHigher = 2.2, HHigher;

            HFunc(alphaLower, HLower);
            HFunc(alphaHigher, HHigher);

            // At the moment this decision is based on a few quick test cases
            // If the root is not near 2.0 then f - f_eq is too small to make the value of alpha matter
            // Chosen to return 2.0 as that is the LBGK case
            if (HLower * HHigher >= 0.0)
              return 2.0;

            return (hemelb::util::NumericalMethods::Brent(&HFunc,
                                                          alphaLower,
                                                          alphaHigher,
                                                          1.0E-3,
                                                          1.0E-10));
          }
          else
          {
            // This is in case previous Alpha was calculated to be zero (does happen occasionally if f_eq - f is small
            prevAlpha = (prevAlpha < 1.8
              ? 2.0
              : prevAlpha);

            // Accuracy is set to 1.0E-3 as this works for difftest.
            return (hemelb::util::NumericalMethods::NewtonRaphson(&HFunc,
                                                                  prevAlpha,
                                                                  1.0E-3,
                                                                  1.0E-10));
          }

        }

      }
    }
  }
}
