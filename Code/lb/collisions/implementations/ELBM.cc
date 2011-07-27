#include "lb/collisions/implementations/ELBM.h"
#include "lb/collisions/implementations/HFunction.h"
#include "util/utilityFunctions.h"

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
        double* ELBM::tau;

        void ELBM::createAlphaArray(const size_t size)
        {
          alpha = new double[size];
          for (size_t i = 0; i < size; i++)
            alpha[i] = 2.0;
        }

        void ELBM::setTau(double* t)
        {
          tau = t;
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

        // WARNING: DOES NOT CALCULATE ALPHA, BECAUSE NON OF THE CURRENT BCS USE IT
        void ELBM::getBoundarySiteValues(const distribn_t* f,
                                         const distribn_t &density,
                                         const distribn_t &v_x,
                                         const distribn_t &v_y,
                                         const distribn_t &v_z,
                                         distribn_t* f_eq,
                                         const site_t index)
        {
          //currentAlphaIndex = index;
          D3Q15::CalculateEntropicFeq(density, v_x, v_y, v_z, f_eq);
          //alpha[index] = getAlpha(f, f_eq, alpha[index]);
        }

        distribn_t ELBM::getOperatorElement(distribn_t &f_i,
                                            distribn_t &f_neq_i,
                                            const LbmParameters* iLbmParams)
        {
          return (alpha[currentAlphaIndex] * iLbmParams->Beta * f_neq_i);
        }

        double ELBM::getAlpha(const distribn_t* lF, const distribn_t* lFEq, double prevAlpha)
        {
          bool big = false;
          double deviation = 0.0;

          for (unsigned int i = 0; i < D3Q15::NUMVECTORS; i++)
          {
            // Papers suggest f_eq - f < 0.001 or (f_eq - f)/f < 0.01 for the point to have approx alpha = 2
            // Accuracy can change depending on stability requirements, because the more NR evaluations it skips
            // the more of the simulation is in the LBGK limit.
            deviation = util::NumericalFunctions::max(fabs( (lFEq[i] - lF[i]) / lF[i]), deviation);
            if (deviation > 1.0E-2)
            {
              big = true;
            }
          }

          if (big)
          {
            HFunction HFunc(lF, lFEq);

            // This is in case previous Alpha was calculated to be zero (does happen occasionally if f_eq - f is small
            prevAlpha = (prevAlpha < 1.8
              ? 2.0
              : prevAlpha);

            // Accuracy is set to 1.0E-3 as this works for difftest.
            return (hemelb::util::NumericalMethods::NewtonRaphson(&HFunc, prevAlpha, 1.0E-3));
          }
          else
          {
            HFunction HFunc(lF, lFEq);

            // The bracket is very large, but it should guarantee that a root is enclosed
            double alphaLower = 2.0 * (*tau), HLower;
            double alphaHigher = 2.0 * (*tau) / deviation, HHigher;

            HFunc(alphaLower, HLower);
            HFunc(alphaHigher, HHigher);

            // The root should be enclosed, but in case it isn't return some default
            // Chosen to return 2.0 as that is the LBGK case
            if (HLower * HHigher >= 0.0)
              return 2.0;

            return (hemelb::util::NumericalMethods::Brent(&HFunc, alphaLower, alphaHigher, 1.0E-3, 1.0E-3));
          }

        }

      }
    }
  }
}
