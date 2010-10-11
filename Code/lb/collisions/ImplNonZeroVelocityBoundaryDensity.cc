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

  f = &f_old[i * 15];

  for (l = 0; l < 15; l++)
  {
    f_neq[l] = f[l];
  }
  boundary_id = (net->net_site_data[i] & BOUNDARY_ID_MASK) >> BOUNDARY_ID_SHIFT;

  *density = mBoundaryDensityArray[boundary_id];

  DensityAndVelocity(f, &dummy_density, v_x, v_y, v_z);
  CalculateFeq(*density, *v_x, *v_y, *v_z, f);

  for (l = 0; l < 15; l++)
  {
    f_neq[l] -= (f_new[f_id[i * 15 + l]] = f[l]);
  }
}

}
}
}
