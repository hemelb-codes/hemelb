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

        double CollisionOperator::getAlpha(const distribn_t* lF,const  distribn_t* lFEq)
        {
          bool small = true;
          distribn_t deviation[D3Q15::NUMVECTORS];

          for (unsigned int i = 0; i < D3Q15::NUMVECTORS; i++)
          {
            if (fabs(/* deviation[i] = */(lFEq[i] - lF[i]) / lF[i]) > 0.01)
            {
              small = false;
              break;
            }
          }

          if (small)
            return 2.0;

          /*
           if (small)
           {
           double a1 = 0.0, a2 = 0.0, a3 = 0.0, a4 = 0.0;
           double deviation2[D3Q15::NUMVECTORS]; //Will be holding deviation[i]^n * (f_eq - f)

           for (unsigned int i = 0; i < D3Q15::NUMVECTORS; i++)
           {
           a1 += (deviation2[i] = deviation[i] * (lFEq[i] - lF[i]));
           a2 += (deviation2[i] *= deviation[i]);
           a3 += (deviation2[i] *= deviation[i]);
           a4 += (deviation2[i] * deviation[i]);
           }

           a1 *= 0.5;
           a2 *= -1.0 / 6.0;
           a3 *= 1.0 / 12.0;
           a4 *= -0.05;

           // Alpha tends to 2 as deviation -> 0. The evaluation of alpha below fails if a1 is 0 so to prevent
           // NaNs and to save some computation we just return 2
           if (a1 < 1.0E-10)
           return 2.0;

           return (2.0 + 4.0 * (4 * a1 * a2 * (a2 + 5.0 * a3) - a1 * a1 * (a2 + 2.0 * a3
           + 4.0 * a4) - 20.0 * a2 * a2 * a2) / (a1 * a1 * a1));
           }
           */

          HFunction HFunc(lF, lFEq);

          return (hemelb::util::NumericalFunctions::NewtonRaphson(&HFunc, 2.0, 1.0E-6));

        }

      }
    }
  }
}

