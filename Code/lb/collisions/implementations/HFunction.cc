#include "lb/collisions/implementations/HFunction.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {
      namespace implementations
      {

        HFunction::HFunction(const distribn_t* lF, const distribn_t* lFEq) :
          mF(lF), mFEq(lFEq)
        {
        }

        void HFunction::operator()(const double alpha, double &H, double &dH)
        {
          double f_alpha[D3Q15::NUMVECTORS];

          CalculateFalphaAndHInternal(alpha, f_alpha, H);

          dH = 0.0;

          for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ++ii)
          {
            dH += (f_alpha[ii] < 0.0
              ? -1.0
              : 1.0) * (mFEq[ii] - mF[ii]) * (1.0 + log(fabs(f_alpha[ii]) / D3Q15::EQMWEIGHTS[ii]));
          }
        }

        void HFunction::operator()(const double alpha, double &H)
        {
          double f_alpha[D3Q15::NUMVECTORS];

          CalculateFalphaAndHInternal(alpha, f_alpha, H);
        }

        double HFunction::eval()
        {
          double lRet = 0.0;

          for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ++ii)
          {
            lRet += h(mF[ii], 1.0 / D3Q15::EQMWEIGHTS[ii]);
          }

          return lRet;
        }

        void HFunction::CalculateFalphaAndHInternal(const double alpha,
                                                    double fAlpha[D3Q15::NUMVECTORS],
                                                    double &H)
        {
          for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ++ii)
          {
            fAlpha[ii] = mF[ii] + alpha * (mFEq[ii] - mF[ii]);
          }

          H = 0.0;

          for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ++ii)
          {
            H += h(fabs(fAlpha[ii]), 1.0 / D3Q15::EQMWEIGHTS[ii]) - h(mF[ii], 1.0
                / D3Q15::EQMWEIGHTS[ii]);
          }
        }

        double HFunction::h(double fi, double wi_1)
        {
          // Assumes distributions are unlikely to be 0.0
          // If they go negative stabilityTester catches it anyway
          return fi * log(fi * wi_1);
        }

      }
    }
  }
}
