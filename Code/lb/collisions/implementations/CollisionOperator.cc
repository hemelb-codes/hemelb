#include "lb/collisions/implementations/CollisionOperator.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {
      namespace implementations
      {

        CollisionOperator::CollisionOperator()
        {

        }

        double CollisionOperator::getAlpha(const distribn_t* lF,
                                           const distribn_t* lFEq,
                                           double prevAlpha)
        {
          bool small = true;

          for (unsigned int i = 0; i < D3Q15::NUMVECTORS; i++)
          {
            // Papers suggest f_eq - f < 0.001 or (f_eq - f)/f < 0.01 for the point to have approx alpha = 2
            // Accuracy can change depending on stability requirements, because the more NR evaluations it skips
            // the more of the simulation is in the LBGK limit.
            if (fabs(lFEq[i] - lF[i]) > 1.0E-3)
            {
              small = false;
              break;
            }
          }

          HFunction HFunc(lF, lFEq);

          if (small)
          {
            double alphaLower = 1.8, HLower;
            double alphaHigher = 2.2, HHigher;

            HFunc(alphaLower, HLower);
            HFunc(alphaHigher, HHigher);

            // At the moment this decision is based on a few quick test cases
            // If the root is not near 2.0 then f - f_eq is too small to make the value of alpha matter
            // Chosen to return 2.0 as that is the LBGK case
            if (HLower * HHigher >= 0.0)
              return 2.0;

            return (hemelb::util::NumericalFunctions::Brent(&HFunc,
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
            return (hemelb::util::NumericalFunctions::NewtonRaphson(&HFunc,
                                                                    prevAlpha,
                                                                    1.0E-3,
                                                                    1.0E-10));
          }

        }

      }
    }
  }
}

