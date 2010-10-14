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

  f = &f_old[i * D3Q15::NUMVECTORS];

  for (l = 0; l < D3Q15::NUMVECTORS; l++)
  {
    f_neq[l] = f[l];
  }
  boundary_id = (net->net_site_data[i] & BOUNDARY_ID_MASK) >> BOUNDARY_ID_SHIFT;

  *density = mBoundaryDensityArray[boundary_id];
  *v_x = *v_y = *v_z = 0.F;

  D3Q15::CalculateFeq(*density, 0.0, 0.0, 0.0, f);

  for(int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
  {
    f_neq[ii] -= f_new[f_id[i * D3Q15::NUMVECTORS]] = f[ii];
  }
}

}
}
}
