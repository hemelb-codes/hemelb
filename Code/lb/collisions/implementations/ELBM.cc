#include "lb/collisions/implementations/ELBM.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {
      namespace implementations
      {

        double ELBM::alpha;

        void ELBM::getSiteValues(const distribn_t* f,
                                 distribn_t &density,
                                 distribn_t &v_x,
                                 distribn_t &v_y,
                                 distribn_t &v_z,
                                 distribn_t* f_eq)
        {
          D3Q15::CalculateEntropicDensityVelocityFEq(f, density, v_x, v_y, v_z, f_eq);
          alpha = getAlpha(f, f_eq);
        }

        void ELBM::getBoundarySiteValues(const distribn_t* f,
                                         const distribn_t &density,
                                         const distribn_t &v_x,
                                         const distribn_t &v_y,
                                         const distribn_t &v_z,
                                         distribn_t* f_eq)
        {
          D3Q15::CalculateEntropicFeq(density, v_x, v_y, v_z, f_eq);
          alpha = getAlpha(f, f_eq);
        }

        // Also updates lFEq_i to be lFNeq_i
        distribn_t ELBM::getOperatorElement(distribn_t &f_i,
                                            distribn_t &f_eq_i,
                                            const LbmParameters* iLbmParams)
        {
          return (alpha * iLbmParams->Beta * (f_eq_i = f_i - f_eq_i));
        }

      }
    }
  }
}
