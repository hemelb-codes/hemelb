#include "lb/collisions/ImplZeroVelocityEquilibrium.h"

namespace hemelb
{
namespace lb
{
namespace collisions
{

void ImplZeroVelocityEquilibrium::DoCollisions(double omega, int i, double *density,
  double *v_x, double *v_y, double *v_z, double f_neq[], Net* net)
{
  double *f;
  double temp;

  int l;

  f = &f_old[i * D3Q15::NUMVECTORS];

  for (l = 0; l < D3Q15::NUMVECTORS; l++)
  {
    f_neq[l] = f[l];
  }
  *v_x = *v_y = *v_z = 0.F;

  *density = 0.;

  for (l = 0; l < D3Q15::NUMVECTORS; l++)
    *density += f[l];

  D3Q15::CalculateFeq(*density, 0.0, 0.0, 0.0, f);

  for(int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
  {
    f_neq[ii] -= (f_new[f_id[ii * D3Q15::NUMVECTORS + ii]]) = f[ii];
  }
}

}
}
}
