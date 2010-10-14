#include "lb/collisions/ImplGuoZhengShi.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      void ImplGuoZhengShi::DoCollisions(double omega, int i, double *density,
        double *v_x, double *v_y, double *v_z, double f_neq[], Net* net)
      {
        // First do a normal collision & streaming step, as if we were mid-fluid.
        // NOTE that we use the version that preserves f_old.
        // NOTE that this handily works out the equilibrium density, v_x, v_y and v_z for us
        double lFEq[15];

        double *f = &f_old[i * D3Q15::NUMVECTORS];

        D3Q15::CalculateDensityVelocityFEq(f, density, v_x, v_y, v_z, lFEq);

        for (int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
        {
          f_new[f_id[i * D3Q15::NUMVECTORS + ii]] = f[ii] + omega * (f_neq[ii] = f[ii]
              - lFEq[ii]);
        }

        // Now fill in the un-streamed-to distributions (those that point away from boundaries).
        unsigned int boundary_config = GetBoundaryConfig(net, i);

        double* cut_dists = & (net->cut_distances[i * (D3Q15::NUMVECTORS)]);

        for (int l = 1; l < D3Q15::NUMVECTORS; l++)
        {
          int lAwayFromWallIndex = D3Q15::INVERSEDIRECTIONS[l];

          if (0 != (boundary_config & (1U << (l - 1))))
          {
            double delta = cut_dists[l - 1];
            double uWall[3];
            double fNeq;

            // Work out uw1 (noting that ub is 0 until we implement moving walls)
            uWall[0] = (1 - 1. / delta) * *v_x;
            uWall[1] = (1 - 1. / delta) * *v_y;
            uWall[2] = (1 - 1. / delta) * *v_z;
            fNeq = f_neq[lAwayFromWallIndex];

            // Interpolate with uw2 if delta < 0.75
            if (delta < 0.75)
            {
              // Only do the extra interpolation if there's gonna be a point there to interpolate from, i.e. there's no boundary
              // in the direction of awayFromWallIndex
              if (0 == (boundary_config & (1U << (lAwayFromWallIndex - 1))))
              {
                // Need some info about the next node away from the wall in this direction...
                int nextIOut = f_id[i * D3Q15::NUMVECTORS + lAwayFromWallIndex]
                    / D3Q15::NUMVECTORS;
                double nextNodeDensity, nextNodeV[3], nextNodeFEq[D3Q15::NUMVECTORS];

                D3Q15::CalculateDensityVelocityFEq(&f_old[nextIOut * D3Q15::NUMVECTORS],
                                                   &nextNodeDensity, &nextNodeV[0],
                                                   &nextNodeV[1], &nextNodeV[2],
                                                   &nextNodeFEq[0]);

                for (int a = 0; a < 3; a++)
                  uWall[a] = delta * uWall[a] + (1. - delta) * (delta - 1.)
                      * nextNodeV[a] / (1. + delta);

                fNeq = delta * fNeq + (1. - delta) * (f_old[nextIOut * D3Q15::NUMVECTORS
                    + lAwayFromWallIndex] - nextNodeFEq[lAwayFromWallIndex]);
              }
              // If there's nothing to extrapolate from we, very lamely, do a 0VE-style operation to fill in the missing velocity.
              else
              {
                for (int a = 0; a < 3; a++)
                  uWall[a] = 0.0;//delta * uWall[a];

                fNeq = 0.0;//delta * fNeq;
              }
            }

            // Use a helper function to calculate the actual value of f_eq in the desired direction at the wall node.
            // Note that we assume that the density is the same as at this node
            double fEqTemp[D3Q15::NUMVECTORS];
            D3Q15::CalculateFeq(*density, uWall[0], uWall[1], uWall[2], fEqTemp);

            // Collide and stream!
            f_new[i * D3Q15::NUMVECTORS + lAwayFromWallIndex]
                = fEqTemp[lAwayFromWallIndex] + (1.0 + omega) * fNeq;
          }
        }
      }

    }
  }
}
