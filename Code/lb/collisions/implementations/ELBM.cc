#include "lb/collisions/implementations/ELBM.h"

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

        void ELBM::createAlphaArray(const size_t size)
        {
          alpha = new double[size];
          for (size_t i = 0; i < size; i++)
            alpha[i] = 2.0;
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

        void ELBM::getBoundarySiteValues(const distribn_t* f,
                                         const distribn_t &density,
                                         const distribn_t &v_x,
                                         const distribn_t &v_y,
                                         const distribn_t &v_z,
                                         distribn_t* f_eq,
                                         const site_t index)
        {
          currentAlphaIndex = index;
          D3Q15::CalculateEntropicFeq(density, v_x, v_y, v_z, f_eq);
          alpha[index] = getAlpha(f, f_eq, alpha[index]);
        }

        // Also updates lFEq_i to be lFNeq_i
        distribn_t ELBM::getOperatorElement(distribn_t &f_i,
                                            distribn_t &f_eq_i,
                                            const LbmParameters* iLbmParams)
        {
          return (alpha[currentAlphaIndex] * iLbmParams->Beta * (f_eq_i = f_i - f_eq_i));
        }

      }
    }
  }
}
