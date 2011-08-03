#ifndef HEMELB_LB_HFUNCTION_H
#define HEMELB_LB_HFUNCTION_H

#include <math.h>
#include <cstdlib>

#include "constants.h"
#include "D3Q15.h"

namespace hemelb
{
  namespace lb
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

      private:
        const distribn_t* mF;
        const distribn_t* mFEq;

        void CalculateFalphaAndHInternal(const double alpha,
                                         double fAlpha[D3Q15::NUMVECTORS],
                                         double &H);
        double h(double fi, double wi_1);
    };

  }
}

#endif /* HEMELB_LB_HFUNCTION_H */
