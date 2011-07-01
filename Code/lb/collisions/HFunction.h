#ifndef HEMELB_LB_COLLISIONS_HFUNCTION_H
#define HEMELB_LB_COLLISIONS_HFUNCTION_H

#include <math.h>
#include "constants.h"
#include <cstdlib>

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {
      /*
       * The HFunction class calculates the H function as the name suggests
       * It is a class as the Newton-Raphson function in util takes in an object
       * with an overloaded () operator.
       */
      class HFunction
      {
        public:
          HFunction(distribn_t* lF, distribn_t* lFEq);

          void operator()(const double alpha, double &H, double &dH);

        private:
          distribn_t* mF;
          distribn_t* mFEq;

          double h(double fi, double wi_1);
      };

    }
  }
}

#endif /* HEMELB_LB_COLLISIONS_HFUNCTION_H */
