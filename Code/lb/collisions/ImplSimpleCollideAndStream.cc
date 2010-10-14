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
        double density_1;
        double *f;
        double v_xx, v_yy, v_zz;

        double lFEq[15];

        f = &f_old[i * D3Q15::NUMVECTORS];

        D3Q15::CalculateDensityVelocityFEq(f, density, v_x, v_y, v_z, lFEq);

        for (int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
        {
          f_new[f_id[i * D3Q15::NUMVECTORS + ii]] = f[ii] += omega * (f_neq[ii] = f[ii]
              - lFEq[ii]);
        }
      }

    }
  }
}

