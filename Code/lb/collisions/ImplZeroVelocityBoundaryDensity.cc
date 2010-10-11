#include "lb/collisions/ImplZeroVelocityBoundaryDensity.h"

namespace hemelb
{
namespace lb
{
namespace collisions
{

ImplZeroVelocityBoundaryDensity::ImplZeroVelocityBoundaryDensity(
                                                               double* iOutletDensityArray)
{
  mBoundaryDensityArray = iOutletDensityArray;
}

// Collision + streaming for fluid lattice sites and adjacent to the outlet and the wall.
void ImplZeroVelocityBoundaryDensity::DoCollisions(double omega, int i, double *density,
                                                  double *v_x, double *v_y, double *v_z,
                                                  double f_neq[], Net* net)
{
  double *f;
  double temp;

  int l;

  unsigned int boundary_id;

  f = &f_old[i * 15];

  for (l = 0; l < 15; l++)
  {
    f_neq[l] = f[l];
  }
  boundary_id = (net->net_site_data[i] & BOUNDARY_ID_MASK) >> BOUNDARY_ID_SHIFT;

  *density = mBoundaryDensityArray[boundary_id];

  *v_x = *v_y = *v_z = 0.F;

  f_neq[0] -= (f_new[f_id[i * 15]] = f[0] = (2.0 / 9.0) * *density);

  temp = (1.0 / 9.0) * *density;

  for (l = 1; l < 7; l++)
    f_neq[l] -= (f_new[f_id[i * 15 + l]] = f[l] = temp);

  temp *= (1.0 / 8.0);

  for (l = 7; l < 15; l++)
    f_neq[l] -= (f_new[f_id[i * 15 + l]] = f[l] = temp);
}

}
}
}
