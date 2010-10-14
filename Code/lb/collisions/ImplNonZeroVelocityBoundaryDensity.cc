#include "lb/collisions/ImplNonZeroVelocityBoundaryDensity.h"

namespace hemelb
{
namespace lb
{
namespace collisions
{

ImplNonZeroVelocityBoundaryDensity::ImplNonZeroVelocityBoundaryDensity(
                                                               double* iOutletDensityArray)
{
  mBoundaryDensityArray = iOutletDensityArray;
}


void ImplNonZeroVelocityBoundaryDensity::DoCollisions(double omega, int i, double *density,
                                                  double *v_x, double *v_y, double *v_z,
                                                  double f_neq[], Net* net)
{
  double *f;
  double dummy_density;

  unsigned int boundary_id, l;

  f = &f_old[i * D3Q15::NUMVECTORS];

  for (l = 0; l < D3Q15::NUMVECTORS; l++)
  {
    f_neq[l] = f[l];
  }
  boundary_id = (net->net_site_data[i] & BOUNDARY_ID_MASK) >> BOUNDARY_ID_SHIFT;

  *density = mBoundaryDensityArray[boundary_id];

  D3Q15::CalculateDensityAndVelocity(f, &dummy_density, v_x, v_y, v_z);
  D3Q15::CalculateFeq(*density, *v_x, *v_y, *v_z, f);

  for (l = 0; l < D3Q15::NUMVECTORS; l++)
  {
    f_neq[l] -= (f_new[f_id[i * D3Q15::NUMVECTORS + l]] = f[l]);
  }
}

}
}
}
