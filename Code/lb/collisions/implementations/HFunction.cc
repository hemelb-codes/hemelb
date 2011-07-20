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
          double f_alpha[15] = { mF[0] + alpha * (mFEq[0] - mF[0]),
                                 mF[1] + alpha * (mFEq[1] - mF[1]),
                                 mF[2] + alpha * (mFEq[2] - mF[2]),
                                 mF[3] + alpha * (mFEq[3] - mF[3]),
                                 mF[4] + alpha * (mFEq[4] - mF[4]),
                                 mF[5] + alpha * (mFEq[5] - mF[5]),
                                 mF[6] + alpha * (mFEq[6] - mF[6]),
                                 mF[7] + alpha * (mFEq[7] - mF[7]),
                                 mF[8] + alpha * (mFEq[8] - mF[8]),
                                 mF[9] + alpha * (mFEq[9] - mF[9]),
                                 mF[10] + alpha * (mFEq[10] - mF[10]),
                                 mF[11] + alpha * (mFEq[11] - mF[11]),
                                 mF[12] + alpha * (mFEq[12] - mF[12]),
                                 mF[13] + alpha * (mFEq[13] - mF[13]),
                                 mF[14] + alpha * (mFEq[14] - mF[14]) };

          H = h(fabs(f_alpha[0]), 4.5) - h(mF[0], 4.5) + h(fabs(f_alpha[1]), 9.0) - h(mF[1], 9.0)
              + h(fabs(f_alpha[2]), 9.0) - h(mF[2], 9.0) + h(fabs(f_alpha[3]), 9.0) - h(mF[3], 9.0)
              + h(fabs(f_alpha[4]), 9.0) - h(mF[4], 9.0) + h(fabs(f_alpha[5]), 9.0) - h(mF[5], 9.0)
              + h(fabs(f_alpha[6]), 9.0) - h(mF[6], 9.0) + h(fabs(f_alpha[7]), 72.0) - h(mF[7],
                                                                                         72.0)
              + h(fabs(f_alpha[8]), 72.0) - h(mF[8], 72.0) + h(fabs(f_alpha[9]), 72.0) - h(mF[9],
                                                                                           72.0)
              + h(fabs(f_alpha[10]), 72.0) - h(mF[10], 72.0) + h(fabs(f_alpha[11]), 72.0)
              - h(mF[11], 72.0) + h(fabs(f_alpha[12]), 72.0) - h(mF[12], 72.0)
              + h(fabs(f_alpha[13]), 72.0) - h(mF[13], 72.0) + h(fabs(f_alpha[14]), 72.0)
              - h(mF[14], 72.0);

          dH = (f_alpha[0] < 0.0
            ? -1.0
            : 1.0) * (mFEq[0] - mF[0]) * (1.0 + log(4.5 * fabs(f_alpha[0]))) + (f_alpha[1] < 0.0
            ? -1.0
            : 1.0) * (mFEq[1] - mF[1]) * (1.0 + log(9.0 * fabs(f_alpha[1]))) + (f_alpha[2] < 0.0
            ? -1.0
            : 1.0) * (mFEq[2] - mF[2]) * (1.0 + log(9.0 * fabs(f_alpha[2]))) + (f_alpha[3] < 0.0
            ? -1.0
            : 1.0) * (mFEq[3] - mF[3]) * (1.0 + log(9.0 * fabs(f_alpha[3]))) + (f_alpha[4] < 0.0
            ? -1.0
            : 1.0) * (mFEq[4] - mF[4]) * (1.0 + log(9.0 * fabs(f_alpha[4]))) + (f_alpha[5] < 0.0
            ? -1.0
            : 1.0) * (mFEq[5] - mF[5]) * (1.0 + log(9.0 * fabs(f_alpha[5]))) + (f_alpha[6] < 0.0
            ? -1.0
            : 1.0) * (mFEq[6] - mF[6]) * (1.0 + log(9.0 * fabs(f_alpha[6]))) + (f_alpha[7] < 0.0
            ? -1.0
            : 1.0) * (mFEq[7] - mF[7]) * (1.0 + log(72.0 * fabs(f_alpha[7]))) + (f_alpha[8] < 0.0
            ? -1.0
            : 1.0) * (mFEq[8] - mF[8]) * (1.0 + log(72.0 * fabs(f_alpha[8]))) + (f_alpha[9] < 0.0
            ? -1.0
            : 1.0) * (mFEq[9] - mF[9]) * (1.0 + log(72.0 * fabs(f_alpha[9]))) + (f_alpha[10] < 0.0
            ? -1.0
            : 1.0) * (mFEq[10] - mF[10]) * (1.0 + log(72.0 * fabs(f_alpha[10]))) + (f_alpha[11]
              < 0.0
            ? -1.0
            : 1.0) * (mFEq[11] - mF[11]) * (1.0 + log(72.0 * fabs(f_alpha[11]))) + (f_alpha[12]
              < 0.0
            ? -1.0
            : 1.0) * (mFEq[12] - mF[12]) * (1.0 + log(72.0 * fabs(f_alpha[12]))) + (f_alpha[13]
              < 0.0
            ? -1.0
            : 1.0) * (mFEq[13] - mF[13]) * (1.0 + log(72.0 * fabs(f_alpha[13]))) + (f_alpha[14]
              < 0.0
            ? -1.0
            : 1.0) * (mFEq[14] - mF[14]) * (1.0 + log(72.0 * fabs(f_alpha[14])));

        }

        void HFunction::operator()(const double alpha, double &H)
        {
          double f_alpha[15] = { mF[0] + alpha * (mFEq[0] - mF[0]),
                                 mF[1] + alpha * (mFEq[1] - mF[1]),
                                 mF[2] + alpha * (mFEq[2] - mF[2]),
                                 mF[3] + alpha * (mFEq[3] - mF[3]),
                                 mF[4] + alpha * (mFEq[4] - mF[4]),
                                 mF[5] + alpha * (mFEq[5] - mF[5]),
                                 mF[6] + alpha * (mFEq[6] - mF[6]),
                                 mF[7] + alpha * (mFEq[7] - mF[7]),
                                 mF[8] + alpha * (mFEq[8] - mF[8]),
                                 mF[9] + alpha * (mFEq[9] - mF[9]),
                                 mF[10] + alpha * (mFEq[10] - mF[10]),
                                 mF[11] + alpha * (mFEq[11] - mF[11]),
                                 mF[12] + alpha * (mFEq[12] - mF[12]),
                                 mF[13] + alpha * (mFEq[13] - mF[13]),
                                 mF[14] + alpha * (mFEq[14] - mF[14]) };

          H = h(fabs(f_alpha[0]), 4.5) - h(mF[0], 4.5) + h(fabs(f_alpha[1]), 9.0) - h(mF[1], 9.0)
              + h(fabs(f_alpha[2]), 9.0) - h(mF[2], 9.0) + h(fabs(f_alpha[3]), 9.0) - h(mF[3], 9.0)
              + h(fabs(f_alpha[4]), 9.0) - h(mF[4], 9.0) + h(fabs(f_alpha[5]), 9.0) - h(mF[5], 9.0)
              + h(fabs(f_alpha[6]), 9.0) - h(mF[6], 9.0) + h(fabs(f_alpha[7]), 72.0) - h(mF[7],
                                                                                         72.0)
              + h(fabs(f_alpha[8]), 72.0) - h(mF[8], 72.0) + h(fabs(f_alpha[9]), 72.0) - h(mF[9],
                                                                                           72.0)
              + h(fabs(f_alpha[10]), 72.0) - h(mF[10], 72.0) + h(fabs(f_alpha[11]), 72.0)
              - h(mF[11], 72.0) + h(fabs(f_alpha[12]), 72.0) - h(mF[12], 72.0)
              + h(fabs(f_alpha[13]), 72.0) - h(mF[13], 72.0) + h(fabs(f_alpha[14]), 72.0)
              - h(mF[14], 72.0);
        }

        double HFunction::eval()
        {
          return h(mF[0], 4.5) + h(mF[1], 9.0) + h(mF[2], 9.0) + h(mF[3], 9.0) + h(mF[4], 9.0)
              + h(mF[5], 9.0) + h(mF[6], 9.0) + h(mF[7], 72.0) + h(mF[8], 72.0) + h(mF[9], 72.0)
              + h(mF[10], 72.0) + h(mF[11], 72.0) + h(mF[12], 72.0) + h(mF[13], 72.0) + h(mF[14],
                                                                                          72.0);
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
}
