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
  double temp1, temp2;

  f = &f_old[i * 15];

  *v_x = f[1] + (f[7] + f[9]) + (f[11] + f[13]);
  *v_y = f[3] + (f[12] + f[14]);
  *v_z = f[5] + f[10];

  *density = f[0] + (f[2] + f[4]) + (f[6] + f[8]) + *v_x + *v_y + *v_z;

  *v_x -= f[2] + (f[8] + f[10]) + (f[12] + f[14]);
  *v_y += (f[7] + f[9]) - (f[4] + (f[8] + f[10]) + (f[11] + f[13]));
  *v_z += f[7] + f[11] + f[14] - ( (f[6] + f[8]) + f[9] + f[12] + f[13]);

  density_1 = 1. / *density;

  v_xx = *v_x * *v_x;
  v_yy = *v_y * *v_y;
  v_zz = *v_z * *v_z;

  f_new[f_id[i * 15 + 0]] = f[0] + omega * (f_neq[0] = f[0] - ( (2.0 / 9.0) * *density
      - (1.0 / 3.0) * ( (v_xx + v_yy + v_zz) * density_1)));

  temp1 = (1.0 / 9.0) * *density - (1.0 / 6.0) * ( (v_xx + v_yy + v_zz) * density_1);

  f_new[f_id[i * 15 + 1]] = f[1] += omega * (f_neq[1] = f[1] - ( (temp1 + (0.5
      * density_1) * v_xx) + ( (1.0 / 3.0) * *v_x))); // (+1, 0, 0)
  f_new[f_id[i * 15 + 2]] = f[2] += omega * (f_neq[2] = f[2] - ( (temp1 + (0.5
      * density_1) * v_xx) - ( (1.0 / 3.0) * *v_x))); // (+1, 0, 0)

  f_new[f_id[i * 15 + 3]] = f[3] += omega * (f_neq[3] = f[3] - ( (temp1 + (0.5
      * density_1) * v_yy) + ( (1.0 / 3.0) * *v_y))); // (0, +1, 0)
  f_new[f_id[i * 15 + 4]] = f[4] += omega * (f_neq[4] = f[4] - ( (temp1 + (0.5
      * density_1) * v_yy) - ( (1.0 / 3.0) * *v_y))); // (0, -1, 0)

  f_new[f_id[i * 15 + 5]] = f[5] += omega * (f_neq[5] = f[5] - ( (temp1 + (0.5
      * density_1) * v_zz) + ( (1.0 / 3.0) * *v_z))); // (0, 0, +1)
  f_new[f_id[i * 15 + 6]] = f[6] += omega * (f_neq[6] = f[6] - ( (temp1 + (0.5
      * density_1) * v_zz) - ( (1.0 / 3.0) * *v_z))); // (0, 0, -1)

  temp1 *= (1.0 / 8.0);

  temp2 = (*v_x + *v_y) + *v_z;

  f_new[f_id[i * 15 + 7]] = f[7] += omega * (f_neq[7] = f[7] - ( (temp1 + ( (1.0 / 16.0)
      * density_1) * temp2 * temp2) + ( (1.0 / 24.0) * temp2))); // (+1, +1, +1)
  f_new[f_id[i * 15 + 8]] = f[8] += omega * (f_neq[8] = f[8] - ( (temp1 + ( (1.0 / 16.0)
      * density_1) * temp2 * temp2) - ( (1.0 / 24.0) * temp2))); // (-1, -1, -1)

  temp2 = (*v_x + *v_y) - *v_z;

  f_new[f_id[i * 15 + 9]] = f[9] += omega * (f_neq[9] = f[9] - ( (temp1 + ( (1.0 / 16.0)
      * density_1) * temp2 * temp2) + ( (1.0 / 24.0) * temp2))); // (+1, +1, -1)
  f_new[f_id[i * 15 + 10]] = f[10] += omega * (f_neq[10] = f[10] - ( (temp1 + ( (1.0
      / 16.0) * density_1) * temp2 * temp2) - ( (1.0 / 24.0) * temp2))); // (-1, -1, +1)

  temp2 = (*v_x - *v_y) + *v_z;

  f_new[f_id[i * 15 + 11]] = f[11] += omega * (f_neq[11] = f[11] - ( (temp1 + ( (1.0
      / 16.0) * density_1) * temp2 * temp2) + ( (1.0 / 24.0) * temp2))); // (+1, -1, +1)
  f_new[f_id[i * 15 + 12]] = f[12] += omega * (f_neq[12] = f[12] - ( (temp1 + ( (1.0
      / 16.0) * density_1) * temp2 * temp2) - ( (1.0 / 24.0) * temp2))); // (-1, +1, -1)

  temp2 = (*v_x - *v_y) - *v_z;

  f_new[f_id[i * 15 + 13]] = f[13] += omega * (f_neq[13] = f[13] - ( (temp1 + ( (1.0
      / 16.0) * density_1) * temp2 * temp2) + ( (1.0 / 24.0) * temp2))); // (+1, -1, -1)
  f_new[f_id[i * 15 + 14]] = f[14] += omega * (f_neq[14] = f[14] - ( (temp1 + ( (1.0
      / 16.0) * density_1) * temp2 * temp2) - ( (1.0 / 24.0) * temp2))); // (-1, +1, +1)
}

}
}
}

