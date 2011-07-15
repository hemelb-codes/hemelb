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
          HFunction(const distribn_t* lF, const distribn_t* lFEq);

          void operator()(const double alpha, double &H, double &dH);
          void operator()(const double alpha, double &H);
          double eval();

          // TODO: Make private once testing finished
        private:
          const distribn_t* mF;
          const distribn_t* mFEq;

          double h(double fi, double wi_1);
      };

    }
  }
}

#endif /* HEMELB_LB_COLLISIONS_HFUNCTION_H */
