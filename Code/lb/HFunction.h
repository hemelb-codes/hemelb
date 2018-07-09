
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_HFUNCTION_H
#define HEMELB_LB_HFUNCTION_H

#include <cmath>
#include <cstdlib>

#include "constants.h"

namespace hemelb
{
  namespace lb
  {
    /*
     * The HFunction class calculates the H function as the name suggests
     * It is a class as the Newton-Raphson function in util takes in an object
     * with an overloaded () operator.
     */
    template<class LatticeType>
    class HFunction
    {
      public:
        HFunction(const distribn_t* lF, const distribn_t* lFEq) :
            mF(lF), mFEq(lFEq)
        {

        }

        void operator()(const double alpha, double &H, double &dH)
        {
          double f_alpha[LatticeType::NUMVECTORS];

          CalculateFalphaAndHInternal(alpha, f_alpha, H);

          dH = 0.0;

          for (unsigned int ii = 0; ii < LatticeType::NUMVECTORS; ++ii)
          {
            dH += (f_alpha[ii] < 0.0 ?
              -1.0 :
              1.0) * (mFEq[ii] - mF[ii]) * (1.0 + std::log(std::fabs(f_alpha[ii]) / LatticeType::EQMWEIGHTS[ii]));
          }
        }

        void operator()(const double alpha, double &H)
        {
          double f_alpha[LatticeType::NUMVECTORS];

          CalculateFalphaAndHInternal(alpha, f_alpha, H);
        }

        double eval()
        {
          double lRet = 0.0;

          for (unsigned int ii = 0; ii < LatticeType::NUMVECTORS; ++ii)
          {
            lRet += h(mF[ii], 1.0 / LatticeType::EQMWEIGHTS[ii]);
          }

          return lRet;
        }

      private:
        const distribn_t* mF;
        const distribn_t* mFEq;

        void CalculateFalphaAndHInternal(const double alpha, double* fAlpha, double &H)
        {
          for (unsigned int ii = 0; ii < LatticeType::NUMVECTORS; ++ii)
          {
            fAlpha[ii] = mF[ii] + alpha * (mFEq[ii] - mF[ii]);
          }

          H = 0.0;

          for (unsigned int ii = 0; ii < LatticeType::NUMVECTORS; ++ii)
          {
            H += h(std::fabs(fAlpha[ii]), 1.0 / LatticeType::EQMWEIGHTS[ii]) - h(mF[ii], 1.0 / LatticeType::EQMWEIGHTS[ii]);
          }
        }

        double h(double fi, double wi_1)
        {
          // Assumes distributions are unlikely to be 0.0
          // If they go negative stabilityTester catches it anyway
          return fi * std::log(fi * wi_1);
        }
    };

  }
}

#endif /* HEMELB_LB_HFUNCTION_H */
