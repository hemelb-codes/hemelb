#include "lb/collisions/ImplFInterpolation.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      void ImplFInterpolation::DoCollisions(double omega, int i, double *density,
        double *v_x, double *v_y, double *v_z, double f_neq[], Net* net)
      {
        double *f = &f_old[i * D3Q15::NUMVECTORS];

        // Temporarily store f_eq in f_neq. Rectified later.
        D3Q15::CalculateDensityVelocityFEq(f, density, v_x, v_y, v_z, f_neq);

        for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
        {
          f_new[f_id[i * D3Q15::NUMVECTORS + ii]] = f[ii] += omega * (f_neq[ii] = f[ii]
              - f_neq[ii]);
        }
      }

      bool ImplFInterpolation::PostStep(double omega, int i, double *density,
        double *v_x, double *v_y, double *v_z, double f_neq[], Net* net)
      {
        unsigned int boundary_config = GetBoundaryConfig(net, i);

        // Handle odd indices, then evens - it's slightly easier to take the odd
        // and even cases separately.
        for (unsigned int l = 1; l < D3Q15::NUMVECTORS; l++)
        {
          // Only do it if there's no boundary in the opposite direction, otherwise, just leave at eqm
          if ( (0 != (boundary_config & (1U << (l - 1)))))// && (0 == (boundary_config & (1U << (l)))))
          {
            double twoQ = 2.0 * net->GetCutDistance(i, l);
            if (twoQ < 1.0)
            {
              f_new[i * D3Q15::NUMVECTORS + D3Q15::INVERSEDIRECTIONS[l]] = f_new[i
                  * D3Q15::NUMVECTORS + l] + twoQ * (f_old[i * D3Q15::NUMVECTORS + l]
                  - f_new[i * D3Q15::NUMVECTORS + l]);
            }
            else
            {
              f_new[i * D3Q15::NUMVECTORS + D3Q15::INVERSEDIRECTIONS[l]] = f_old[i
                  * D3Q15::NUMVECTORS + D3Q15::INVERSEDIRECTIONS[l]] + (1. / twoQ)
                  * (f_old[i * D3Q15::NUMVECTORS + l] - f_old[i * 15
                      + D3Q15::INVERSEDIRECTIONS[l]]);
            }
          }
        }

        // Return true because we've done a post step.
        return true;
      }

    }
  }
}
