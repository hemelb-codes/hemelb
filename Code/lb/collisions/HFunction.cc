#include "lb/collisions/HFunction.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      HFunction::HFunction(const distribn_t* lF, const distribn_t* lFEq) :
        mF(lF), mFEq(lFEq)
      {
      }

      void HFunction::operator()(const double alpha, double &H, double &dH)
      {
        double f_alpha = mF[0] + alpha * (mFEq[0] - mF[0]);
        H = h(fabs(f_alpha), 9.0 / 2.0) - h(mF[0], 9.0 / 2.0);
        dH = (f_alpha < 0.0
          ? -1.0
          : 1.0) * (mFEq[0] - mF[0]) * (1.0 + log( (9.0 / 2.0) * fabs(f_alpha)));

        for (int i = 1; i < 7; i++)
        {
          f_alpha = mF[i] + alpha * (mFEq[i] - mF[i]);
          H += h(fabs(f_alpha), 9.0) - h(mF[i], 9.0);
          dH += (f_alpha < 0.0
            ? -1.0
            : 1.0) * (mFEq[i] - mF[i]) * (1.0 + log(9.0 * fabs(f_alpha)));
        }

        for (int i = 7; i < 15; i++)
        {
          f_alpha = mF[i] + alpha * (mFEq[i] - mF[i]);
          H += h(fabs(f_alpha), 72.0) - h(mF[i], 72.0);
          dH += (f_alpha < 0.0
            ? -1.0
            : 1.0) * (mFEq[i] - mF[i]) * (1.0 + log(72.0 * fabs(f_alpha)));
        }
      }

      void HFunction::operator()(const double alpha, double &H)
      {
        double f_alpha = mF[0] + alpha * (mFEq[0] - mF[0]);
        H = h(fabs(f_alpha), 9.0 / 2.0) - h(mF[0], 9.0 / 2.0);

        for (int i = 1; i < 7; i++)
        {
          f_alpha = mF[i] + alpha * (mFEq[i] - mF[i]);
          H += h(fabs(f_alpha), 9.0) - h(mF[i], 9.0);
        }

        for (int i = 7; i < 15; i++)
        {
          f_alpha = mF[i] + alpha * (mFEq[i] - mF[i]);
          H += h(fabs(f_alpha), 72.0) - h(mF[i], 72.0);
        }
      }

      double HFunction::evaluate()
      {
        double H = h(mF[0], 9.0 / 2.0);

        for (int i = 1; i < 7; i++)
        {
          H += h(mF[i], 9.0);
        }

        for (int i = 7; i < 15; i++)
        {
          H += h(mF[i], 72.0);
        }

        return H;
      }

      double HFunction::h(double fi, double wi_1)
      {
        return (fi > 1.0E-10
          ? (fi * log(fi * wi_1))
          : 0.0);
      }

    }
  }
}
