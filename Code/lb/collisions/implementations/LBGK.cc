#include "lb/collisions/implementations/LBGK.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {
      namespace implementations
      {

        void LBGK::getSiteValues(const distribn_t* f,
                                 distribn_t &density,
                                 distribn_t &v_x,
                                 distribn_t &v_y,
                                 distribn_t &v_z,
                                 distribn_t* f_eq)
        {
          D3Q15::CalculateDensityVelocityFEq(f, density, v_x, v_y, v_z, f_eq);
        }

        void LBGK::getBoundarySiteValues(const distribn_t* f,
                                         const distribn_t &density,
                                         const distribn_t &v_x,
                                         const distribn_t &v_y,
                                         const distribn_t &v_z,
                                         distribn_t* f_eq)
        {
          D3Q15::CalculateFeq(density, v_x, v_y, v_z, f_eq);
        }

        // Also updates f_eq to be f_i
        distribn_t LBGK::getOperatorElement(distribn_t &f_i,
                                            distribn_t &f_eq_i,
                                            const LbmParameters* iLbmParams)
        {
          return (iLbmParams->Omega * (f_eq_i = f_i - f_eq_i));
        }

      }
    }
  }
}
