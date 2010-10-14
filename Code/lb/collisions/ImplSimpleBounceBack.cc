#include "lb/collisions/ImplSimpleBounceBack.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      void ImplSimpleBounceBack::DoCollisions(double iOmega, int iIndex,
        double* oDensity, double *oVx, double *oVy, double *oVz, double *oFNeq, Net* iNet)
      {
        double *f = &f_old[iIndex * D3Q15::NUMVECTORS];

        // Temporarily store f_eq in the oFNeq array
        D3Q15::CalculateDensityVelocityFEq(f, oDensity, oVx, oVy, oVz, oFNeq);

        unsigned int boundary_config = GetBoundaryConfig(iNet, iIndex);

        for (int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
        {
          // The actual bounce-back lines, including streaming and collision. Basically swap the non-equilibrium components of f in each of the opposing pairs of directions.
          int lStreamTo = (0 == (boundary_config & (1U << (ii - 1))))
            ? f_id[iIndex * D3Q15::NUMVECTORS + ii]
            : ( (iIndex * D3Q15::NUMVECTORS) + D3Q15::INVERSEDIRECTIONS[ii]);

          // Remember, oFNeq currently hold the equilibrium distribution. We
          // simultaneously use this and correct it, here.
          f_new[lStreamTo] = f[ii] + iOmega * (oFNeq[ii] = f[ii] - oFNeq[ii]);
        }
      }

    }
  }
}
