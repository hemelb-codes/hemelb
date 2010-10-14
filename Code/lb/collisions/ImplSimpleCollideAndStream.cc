#include "lb/collisions/ImplSimpleCollideAndStream.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      void ImplSimpleCollideAndStream::DoCollisions(double omega, int i, double *density,
        double *v_x, double *v_y, double *v_z, double f_neq[], Net* net)
      {
        double *f = &f_old[i * D3Q15::NUMVECTORS];

        // Temporarily store f_eq in f_neq (rectified in next statement)
        D3Q15::CalculateDensityVelocityFEq(f, density, v_x, v_y, v_z, f_neq);

        for (int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
        {
          f_new[f_id[i * D3Q15::NUMVECTORS + ii]] = f[ii] += omega * (f_neq[ii] = f[ii]
              - f_neq[ii]);
        }
      }

    }
  }
}

