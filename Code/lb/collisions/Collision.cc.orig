#include "lb/collisions/Collision.h"

namespace hemelb
{
namespace lb
{
namespace collisions
{

// TODO: Duplicate of the one in lb.cc. Don't really want both. Probably should be a property of a LatticeGeometry object
void Collision::DensityAndVelocity(double f[], double *density, double *v_x, double *v_y,
                                   double *v_z)
{
  *v_x = f[1] + (f[7] + f[9]) + (f[11] + f[13]);
  *v_y = f[3] + (f[12] + f[14]);
  *v_z = f[5] + f[10];

  *density = f[0] + (f[2] + f[4]) + (f[6] + f[8]) + *v_x + *v_y + *v_z;

  *v_x -= (f[2] + f[8] + f[10] + (f[12] + f[14]));
  *v_y += (f[7] + f[9]) - ( (f[4] + f[8] + f[10] + (f[11] + f[13])));
  *v_z += f[7] + f[11] + f[14] - ( ( (f[6] + f[8]) + f[9] + f[12] + f[13]));
}

// TODO: Duplicate of one in lb.cc. Should probably be part of some LatticeGeometry object.
void Collision::CalculateFeq(double density, double v_x, double v_y, double v_z,
                             double f_eq[])
{
  double density_1;
  double v_xx, v_yy, v_zz;
  double temp1, temp2;

  density_1 = 1. / density;

  v_xx = v_x * v_x;
  v_yy = v_y * v_y;
  v_zz = v_z * v_z;

  f_eq[0] = (2.0 / 9.0) * density - (1.0 / 3.0) * ( (v_xx + v_yy + v_zz) * density_1);

  temp1 = (1.0 / 9.0) * density - (1.0 / 6.0) * ( (v_xx + v_yy + v_zz) * density_1);

  f_eq[1] = (temp1 + (0.5 * density_1) * v_xx) + ( (1.0 / 3.0) * v_x); // (+1, 0, 0)
  f_eq[2] = (temp1 + (0.5 * density_1) * v_xx) - ( (1.0 / 3.0) * v_x); // (+1, 0, 0)

  f_eq[3] = (temp1 + (0.5 * density_1) * v_yy) + ( (1.0 / 3.0) * v_y); // (0, +1, 0)
  f_eq[4] = (temp1 + (0.5 * density_1) * v_yy) - ( (1.0 / 3.0) * v_y); // (0, -1, 0)

  f_eq[5] = (temp1 + (0.5 * density_1) * v_zz) + ( (1.0 / 3.0) * v_z); // (0, 0, +1)
  f_eq[6] = (temp1 + (0.5 * density_1) * v_zz) - ( (1.0 / 3.0) * v_z); // (0, 0, -1)

  temp1 *= (1.0 / 8.0);

  temp2 = (v_x + v_y) + v_z;

  f_eq[7] = (temp1 + ( (1.0 / 16.0) * density_1) * temp2 * temp2) + ( (1.0 / 24.0)
      * temp2); // (+1, +1, +1)
  f_eq[8] = (temp1 + ( (1.0 / 16.0) * density_1) * temp2 * temp2) - ( (1.0 / 24.0)
      * temp2); // (-1, -1, -1)

  temp2 = (v_x + v_y) - v_z;

  f_eq[9] = (temp1 + ( (1.0 / 16.0) * density_1) * temp2 * temp2) + ( (1.0 / 24.0)
      * temp2); // (+1, +1, -1)
  f_eq[10] = (temp1 + ( (1.0 / 16.0) * density_1) * temp2 * temp2) - ( (1.0 / 24.0)
      * temp2); // (-1, -1, +1)

  temp2 = (v_x - v_y) + v_z;

  f_eq[11] = (temp1 + ( (1.0 / 16.0) * density_1) * temp2 * temp2) + ( (1.0 / 24.0)
      * temp2); // (+1, -1, +1)
  f_eq[12] = (temp1 + ( (1.0 / 16.0) * density_1) * temp2 * temp2) - ( (1.0 / 24.0)
      * temp2); // (-1, +1, -1)

  temp2 = (v_x - v_y) - v_z;

  f_eq[13] = (temp1 + ( (1.0 / 16.0) * density_1) * temp2 * temp2) + ( (1.0 / 24.0)
      * temp2); // (+1, -1, -1)
  f_eq[14] = (temp1 + ( (1.0 / 16.0) * density_1) * temp2 * temp2) - ( (1.0 / 24.0)
      * temp2); // (-1, +1, +1)
}

}
}
}
