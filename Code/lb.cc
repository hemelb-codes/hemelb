// In this file, the functions useful to calculate the equilibrium distribution
// function, momentums, the effective von Mises stress and the boundary conditions
// are reported

#include "config.h"


void (*lbmInnerCollision[COLLISION_TYPES]) (double omega, int i, double *density, double *v_x, double *v_y, double *v_z, double f_neq[]);
void (*lbmInnerCollisionSIMD[COLLISION_TYPES]) (double omega, int i, double density[], double v_x[], double v_y[], double v_z[], double f_neq[]);

void (*lbmInterCollision[COLLISION_TYPES]) (double omega, int i, double *density, double *v_x, double *v_y, double *v_z, double f_neq[]);
void (*lbmInterCollisionSIMD[COLLISION_TYPES]) (double omega, int i, double density[], double v_x[], double v_y[], double v_z[], double f_neq[]);


void lbmFeq (double f[], double *density, double *v_x, double *v_y, double *v_z, double f_eq[])
{
  double v_xx, v_yy, v_zz;
  double temp1, temp2;
  
  
  *v_x = f[1] + (f[7] + f[9]) + (f[11] + f[13]);
  *v_y = f[3] + (f[12] + f[14]);
  *v_z = f[5] + f[10];
  
  *density = f[0] + (f[2] + f[4]) + (f[6] + f[8]) + *v_x + *v_y + *v_z;
  
  *v_x -= f[2] + (f[8] + f[10]) + (f[12] + f[14]);
  *v_y += (f[7] + f[9]) - (f[4] + (f[8] + f[10]) + (f[11] + f[13]));
  *v_z += f[7] + f[11] + f[14] - ((f[6] + f[8]) + f[9] + f[12] + f[13]);
  
  v_xx = *v_x * *v_x;
  v_yy = *v_y * *v_y;
  v_zz = *v_z * *v_z;
  
  temp1 = (1.0 / 8.0) * *density;
  
  temp2 = v_xx + v_yy + v_zz;
  
  f_eq[0] = temp1 - (1.0 / 3.0) * temp2;
  
  temp1 -= (1.0 / 6.0) * temp2;
  
  f_eq[1] = (temp1 + 0.5 * v_xx) + ((1.0 / 3.0) * *v_x);   // (+1, 0, 0)
  f_eq[2] = (temp1 + 0.5 * v_xx) - ((1.0 / 3.0) * *v_x);   // (+1, 0, 0)
  
  f_eq[3] = (temp1 + 0.5 * v_yy) + ((1.0 / 3.0) * *v_y);   // (0, +1, 0)
  f_eq[4] = (temp1 + 0.5 * v_yy) - ((1.0 / 3.0) * *v_y);   // (0, -1, 0)
  
  f_eq[5] = (temp1 + 0.5 * v_zz) + ((1.0 / 3.0) * *v_z);   // (0, 0, +1)
  f_eq[6] = (temp1 + 0.5 * v_zz) - ((1.0 / 3.0) * *v_z);   // (0, 0, -1)
  
  temp1 *= (1.0 / 8.0);
  
  temp2 = (*v_x + *v_y) + *v_z;
  
  f_eq[ 7] = (temp1 + (1.0 / 16.0) * temp2 * temp2) + ((1.0 / 24.0) * temp2);   // (+1, +1, +1)
  f_eq[ 8] = (temp1 + (1.0 / 16.0) * temp2 * temp2) - ((1.0 / 24.0) * temp2);   // (-1, -1, -1)
  
  temp2 = (*v_x + *v_y) - *v_z;
  
  f_eq[ 9] = (temp1 + (1.0 / 16.0) * temp2 * temp2) + ((1.0 / 24.0) * temp2);   // (+1, +1, -1)
  f_eq[10] = (temp1 + (1.0 / 16.0) * temp2 * temp2) - ((1.0 / 24.0) * temp2);   // (-1, -1, +1)
  
  temp2 = (*v_x - *v_y) + *v_z;
  
  f_eq[11] = (temp1 + (1.0 / 16.0) * temp2 * temp2) + ((1.0 / 24.0) * temp2);   // (+1, -1, +1)
  f_eq[12] = (temp1 + (1.0 / 16.0) * temp2 * temp2) - ((1.0 / 24.0) * temp2);   // (-1, +1, -1)
  
  temp2 = (*v_x - *v_y) - *v_z;
  
  f_eq[13] = (temp1 + (1.0 / 16.0) * temp2 * temp2) + ((1.0 / 24.0) * temp2);   // (+1, -1, -1)
  f_eq[14] = (temp1 + (1.0 / 16.0) * temp2 * temp2) - ((1.0 / 24.0) * temp2);   // (-1, +1, +1)
}


void lbmFeq (double density, double v_x, double v_y, double v_z, double f_eq[])
{
  double v_xx, v_yy, v_zz;
  double temp1, temp2;
  
  
  v_xx = v_x * v_x;
  v_yy = v_y * v_y;
  v_zz = v_z * v_z;
  
  temp1 = (1.0 / 8.0) * density;
  
  temp2 = v_xx + v_yy + v_zz;
  
  f_eq[0] = temp1 - (1.0 / 3.0) * temp2;
  
  temp1 -= (1.0 / 6.0) * temp2;
  
  f_eq[1] = (temp1 + 0.5 * v_xx) + ((1.0 / 3.0) * v_x);   // (+1, 0, 0)
  f_eq[2] = (temp1 + 0.5 * v_xx) - ((1.0 / 3.0) * v_x);   // (+1, 0, 0)
  
  f_eq[3] = (temp1 + 0.5 * v_yy) + ((1.0 / 3.0) * v_y);   // (0, +1, 0)
  f_eq[4] = (temp1 + 0.5 * v_yy) - ((1.0 / 3.0) * v_y);   // (0, -1, 0)
  
  f_eq[5] = (temp1 + 0.5 * v_zz) + ((1.0 / 3.0) * v_z);   // (0, 0, +1)
  f_eq[6] = (temp1 + 0.5 * v_zz) - ((1.0 / 3.0) * v_z);   // (0, 0, -1)
  
  temp1 *= (1.0 / 8.0);
  
  temp2 = (v_x + v_y) + v_z;
  
  f_eq[ 7] = (temp1 + (1.0 / 16.0) * temp2 * temp2) + ((1.0 / 24.0) * temp2);   // (+1, +1, +1)
  f_eq[ 8] = (temp1 + (1.0 / 16.0) * temp2 * temp2) - ((1.0 / 24.0) * temp2);   // (-1, -1, -1)
  
  temp2 = (v_x + v_y) - v_z;				     
  
  f_eq[ 9] = (temp1 + (1.0 / 16.0) * temp2 * temp2) + ((1.0 / 24.0) * temp2);   // (+1, +1, -1)
  f_eq[10] = (temp1 + (1.0 / 16.0) * temp2 * temp2) - ((1.0 / 24.0) * temp2);   // (-1, -1, +1)
  
  temp2 = (v_x - v_y) + v_z;				     
  
  f_eq[11] = (temp1 + (1.0 / 16.0) * temp2 * temp2) + ((1.0 / 24.0) * temp2);   // (+1, -1, +1)
  f_eq[12] = (temp1 + (1.0 / 16.0) * temp2 * temp2) - ((1.0 / 24.0) * temp2);   // (-1, +1, -1)
  
  temp2 = (v_x - v_y) - v_z;				     
  						     
  f_eq[13] = (temp1 + (1.0 / 16.0) * temp2 * temp2) + ((1.0 / 24.0) * temp2);   // (+1, -1, -1)
  f_eq[14] = (temp1 + (1.0 / 16.0) * temp2 * temp2) - ((1.0 / 24.0) * temp2);   // (-1, +1, +1)
}


void lbmInnerCollision0 (double omega, int i,
			 double *density, double *v_x, double *v_y, double *v_z,
			 double f_neq[])
{
  double *f;
  double v_xx, v_yy, v_zz;
  double temp1, temp2;
  
  
#ifndef TD
  f = &f_old[ i*15 ];
#else
  f = &f_old[ i*30+is_current*15 ];
#endif
  
  *v_x = f[1] + (f[7] + f[9]) + (f[11] + f[13]);
  *v_y = f[3] + (f[12] + f[14]);
  *v_z = f[5] + f[10];
  
  *density = f[0] + (f[2] + f[4]) + (f[6] + f[ 8 ]) + *v_x + *v_y + *v_z;
  
  *v_x -= f[2] + (f[8] + f[10]) + (f[12] + f[14]);
  *v_y += (f[7] + f[9]) - (f[4] + (f[8] + f[10]) + (f[11] + f[13]));
  *v_z += f[7] + f[11] + f[14] - ((f[6] + f[8]) + f[9] + f[12] + f[13]);

  v_xx = *v_x * *v_x;
  v_yy = *v_y * *v_y;
  v_zz = *v_z * *v_z;
  
  temp1 = (1.0 / 8.0) * *density;
  
  temp2 = v_xx + v_yy + v_zz;
  
#ifndef TD
  f_new[ f_id[i*15+0]               ] = f[0] + omega * (f_neq[0] = f[0] - (temp1 - (1.0 / 3.0) * temp2));
#else
  f_new[ f_id[i*15+0]+is_current*15 ] = f[0] + omega * (f_neq[0] = f[0] - (temp1 - (1.0 / 3.0) * temp2));
#endif
  temp1 -= (1.0 / 6.0) * temp2;
  
#ifndef TD
  f_new[ f_id[i*15+1]               ] = f[1] + omega * (f_neq[1] = f[1] - ((temp1 + 0.5 * v_xx) + ((1.0 / 3.0) * *v_x)));   // (+1, 0, 0)
  f_new[ f_id[i*15+2]               ] = f[2] + omega * (f_neq[2] = f[2] - ((temp1 + 0.5 * v_xx) - ((1.0 / 3.0) * *v_x)));   // (+1, 0, 0)
  
  f_new[ f_id[i*15+3]               ] = f[3] + omega * (f_neq[3] = f[3] - ((temp1 + 0.5 * v_yy) + ((1.0 / 3.0) * *v_y)));   // (0, +1, 0)
  f_new[ f_id[i*15+4]               ] = f[4] + omega * (f_neq[4] = f[4] - ((temp1 + 0.5 * v_yy) - ((1.0 / 3.0) * *v_y)));   // (0, -1, 0)
  
  f_new[ f_id[i*15+5]               ] = f[5] + omega * (f_neq[5] = f[5] - ((temp1 + 0.5 * v_zz) + ((1.0 / 3.0) * *v_z)));   // (0, 0, +1)
  f_new[ f_id[i*15+6]               ] = f[6] + omega * (f_neq[6] = f[6] - ((temp1 + 0.5 * v_zz) - ((1.0 / 3.0) * *v_z)));   // (0, 0, -1)
#else
  f_new[ f_id[i*15+1]+is_current*15 ] = f[1] + omega * (f_neq[1] = f[1] - ((temp1 + 0.5 * v_xx) + ((1.0 / 3.0) * *v_x)));   // (+1, 0, 0)
  f_new[ f_id[i*15+2]+is_current*15 ] = f[2] + omega * (f_neq[2] = f[2] - ((temp1 + 0.5 * v_xx) - ((1.0 / 3.0) * *v_x)));   // (+1, 0, 0)
  
  f_new[ f_id[i*15+3]+is_current*15 ] = f[3] + omega * (f_neq[3] = f[3] - ((temp1 + 0.5 * v_yy) + ((1.0 / 3.0) * *v_y)));   // (0, +1, 0)
  f_new[ f_id[i*15+4]+is_current*15 ] = f[4] + omega * (f_neq[4] = f[4] - ((temp1 + 0.5 * v_yy) - ((1.0 / 3.0) * *v_y)));   // (0, -1, 0)
  
  f_new[ f_id[i*15+5]+is_current*15 ] = f[5] + omega * (f_neq[5] = f[5] - ((temp1 + 0.5 * v_zz) + ((1.0 / 3.0) * *v_z)));   // (0, 0, +1)
  f_new[ f_id[i*15+6]+is_current*15 ] = f[6] + omega * (f_neq[6] = f[6] - ((temp1 + 0.5 * v_zz) - ((1.0 / 3.0) * *v_z)));   // (0, 0, -1)
#endif
  temp1 *= (1.0 / 8.0);
  
  temp2 = (*v_x + *v_y) + *v_z;
  
#ifndef TD
  f_new[ f_id[i*15+7]               ] = f[7] + omega * (f_neq[7] = f[7] - ((temp1 + (1.0 / 16.0) * temp2 * temp2) + ((1.0 / 24.0) * temp2)));   // (+1, +1, +1)
  f_new[ f_id[i*15+8]               ] = f[8] + omega * (f_neq[8] = f[8] - ((temp1 + (1.0 / 16.0) * temp2 * temp2) - ((1.0 / 24.0) * temp2)));   // (-1, -1, -1)
#else
  f_new[ f_id[i*15+7]+is_current*15 ] = f[7] + omega * (f_neq[7] = f[7] - ((temp1 + (1.0 / 16.0) * temp2 * temp2) + ((1.0 / 24.0) * temp2)));   // (+1, +1, +1)
  f_new[ f_id[i*15+8]+is_current*15 ] = f[8] + omega * (f_neq[8] = f[8] - ((temp1 + (1.0 / 16.0) * temp2 * temp2) - ((1.0 / 24.0) * temp2)));   // (-1, -1, -1)
#endif
  
  temp2 = (*v_x + *v_y) - *v_z;
  
#ifndef TD
  f_new[ f_id[i*15+ 9]               ] = f[ 9] + omega * (f_neq[ 9] = f[ 9] - ((temp1 + (1.0 / 16.0) * temp2 * temp2) + ((1.0 / 24.0) * temp2)));   // (+1, +1, -1)
  f_new[ f_id[i*15+10]               ] = f[10] + omega * (f_neq[10] = f[10] - ((temp1 + (1.0 / 16.0) * temp2 * temp2) - ((1.0 / 24.0) * temp2)));   // (-1, -1, +1)
#else
  f_new[ f_id[i*15+ 9]+is_current*15 ] = f[ 9] + omega * (f_neq[ 9] = f[ 9] - ((temp1 + (1.0 / 16.0) * temp2 * temp2) + ((1.0 / 24.0) * temp2)));   // (+1, +1, -1)
  f_new[ f_id[i*15+10]+is_current*15 ] = f[10] + omega * (f_neq[10] = f[10] - ((temp1 + (1.0 / 16.0) * temp2 * temp2) - ((1.0 / 24.0) * temp2)));   // (-1, -1, +1)
#endif
  
  temp2 = (*v_x - *v_y) + *v_z;
  
#ifndef TD
  f_new[ f_id[i*15+11]               ] = f[11] + omega * (f_neq[11] = f[11] - ((temp1 + (1.0 / 16.0) * temp2 * temp2) + ((1.0 / 24.0) * temp2)));   // (+1, -1, +1)
  f_new[ f_id[i*15+12]               ] = f[12] + omega * (f_neq[12] = f[12] - ((temp1 + (1.0 / 16.0) * temp2 * temp2) - ((1.0 / 24.0) * temp2)));   // (-1, +1, -1)
  #else
  f_new[ f_id[i*15+11]+is_current*15 ] = f[11] + omega * (f_neq[11] = f[11] - ((temp1 + (1.0 / 16.0) * temp2 * temp2) + ((1.0 / 24.0) * temp2)));   // (+1, -1, +1)
  f_new[ f_id[i*15+12]+is_current*15 ] = f[12] + omega * (f_neq[12] = f[12] - ((temp1 + (1.0 / 16.0) * temp2 * temp2) - ((1.0 / 24.0) * temp2)));   // (-1, +1, -1)
#endif
  
  temp2 = (*v_x - *v_y) - *v_z;	 
  
#ifndef TD
  f_new[ f_id[i*15+13]               ] = f[13] + omega * (f_neq[13] = f[13] - ((temp1 + (1.0 / 16.0) * temp2 * temp2) + ((1.0 / 24.0) * temp2)));   // (+1, -1, -1)
  f_new[ f_id[i*15+14]               ] = f[14] + omega * (f_neq[14] = f[14] - ((temp1 + (1.0 / 16.0) * temp2 * temp2) - ((1.0 / 24.0) * temp2)));   // (-1, +1, +1)
#else
  f_new[ f_id[i*15+13]+is_current*15 ] = f[13] + omega * (f_neq[13] = f[13] - ((temp1 + (1.0 / 16.0) * temp2 * temp2) + ((1.0 / 24.0) * temp2)));   // (+1, -1, -1)
  f_new[ f_id[i*15+14]+is_current*15 ] = f[14] + omega * (f_neq[14] = f[14] - ((temp1 + (1.0 / 16.0) * temp2 * temp2) - ((1.0 / 24.0) * temp2)));   // (-1, +1, +1)
#endif
  
#ifndef BENCH
  for (int l = 0; l < 15; l++)
    {
      if (f[l] < 0.) is_unstable = 1;
    }
#endif
}


void lbmInterCollision0 (double omega, int i,
			 double *density, double *v_x, double *v_y, double *v_z,
			 double f_neq[])
{
  double *f;
  double v_xx, v_yy, v_zz;
  double temp1, temp2;
  
  
#ifndef TD
  f = &f_old[ i*15 ];
#else
  f = &f_old[ i*30+is_current*15 ];
#endif
  
  *v_x = f[1] + (f[7] + f[9]) + (f[11] + f[13]);
  *v_y = f[3] + (f[12] + f[14]);
  *v_z = f[5] + f[10];
  
  *density = f[0] + (f[2] + f[4]) + (f[6] + f[ 8 ]) + *v_x + *v_y + *v_z;
  
  *v_x -= f[2] + (f[8] + f[10]) + (f[12] + f[14]);
  *v_y += (f[7] + f[9]) - (f[4] + (f[8] + f[10]) + (f[11] + f[13]));
  *v_z += f[7] + f[11] + f[14] - ((f[6] + f[8]) + f[9] + f[12] + f[13]);

  v_xx = *v_x * *v_x;
  v_yy = *v_y * *v_y;
  v_zz = *v_z * *v_z;
  
  temp1 = (1.0 / 8.0) * *density;
  
  temp2 = v_xx + v_yy + v_zz;
  
#ifndef TD
  f_new[ f_id[i*15+0]               ] = f[0] + omega * (f_neq[0] = f[0] - (temp1 - (1.0 / 3.0) * temp2));
#else
  f_new[ f_id[i*15+0]+is_current*15 ] = f[0] + omega * (f_neq[0] = f[0] - (temp1 - (1.0 / 3.0) * temp2));
#endif
  temp1 -= (1.0 / 6.0) * temp2;
  
#ifndef TD
  f_new[ f_id[i*15+1]               ] = f[1] += omega * (f_neq[1] = f[1] - ((temp1 + 0.5 * v_xx) + ((1.0 / 3.0) * *v_x)));   // (+1, 0, 0)
  f_new[ f_id[i*15+2]               ] = f[2] += omega * (f_neq[2] = f[2] - ((temp1 + 0.5 * v_xx) - ((1.0 / 3.0) * *v_x)));   // (+1, 0, 0)
  
  f_new[ f_id[i*15+3]               ] = f[3] += omega * (f_neq[3] = f[3] - ((temp1 + 0.5 * v_yy) + ((1.0 / 3.0) * *v_y)));   // (0, +1, 0)
  f_new[ f_id[i*15+4]               ] = f[4] += omega * (f_neq[4] = f[4] - ((temp1 + 0.5 * v_yy) - ((1.0 / 3.0) * *v_y)));   // (0, -1, 0)
  
  f_new[ f_id[i*15+5]               ] = f[5] += omega * (f_neq[5] = f[5] - ((temp1 + 0.5 * v_zz) + ((1.0 / 3.0) * *v_z)));   // (0, 0, +1)
  f_new[ f_id[i*15+6]               ] = f[6] += omega * (f_neq[6] = f[6] - ((temp1 + 0.5 * v_zz) - ((1.0 / 3.0) * *v_z)));   // (0, 0, -1)
#else
  f_new[ f_id[i*15+1]+is_current*15 ] = f[1] += omega * (f_neq[1] = f[1] - ((temp1 + 0.5 * v_xx) + ((1.0 / 3.0) * *v_x)));   // (+1, 0, 0)
  f_new[ f_id[i*15+2]+is_current*15 ] = f[2] += omega * (f_neq[2] = f[2] - ((temp1 + 0.5 * v_xx) - ((1.0 / 3.0) * *v_x)));   // (+1, 0, 0)
  
  f_new[ f_id[i*15+3]+is_current*15 ] = f[3] += omega * (f_neq[3] = f[3] - ((temp1 + 0.5 * v_yy) + ((1.0 / 3.0) * *v_y)));   // (0, +1, 0)
  f_new[ f_id[i*15+4]+is_current*15 ] = f[4] += omega * (f_neq[4] = f[4] - ((temp1 + 0.5 * v_yy) - ((1.0 / 3.0) * *v_y)));   // (0, -1, 0)
  
  f_new[ f_id[i*15+5]+is_current*15 ] = f[5] += omega * (f_neq[5] = f[5] - ((temp1 + 0.5 * v_zz) + ((1.0 / 3.0) * *v_z)));   // (0, 0, +1)
  f_new[ f_id[i*15+6]+is_current*15 ] = f[6] += omega * (f_neq[6] = f[6] - ((temp1 + 0.5 * v_zz) - ((1.0 / 3.0) * *v_z)));   // (0, 0, -1)
#endif
  temp1 *= (1.0 / 8.0);
  
  temp2 = (*v_x + *v_y) + *v_z;
  
#ifndef TD
  f_new[ f_id[i*15+7]               ] = f[7] += omega * (f_neq[7] = f[7] - ((temp1 + (1.0 / 16.0) * temp2 * temp2) + ((1.0 / 24.0) * temp2)));   // (+1, +1, +1)
  f_new[ f_id[i*15+8]               ] = f[8] += omega * (f_neq[8] = f[8] - ((temp1 + (1.0 / 16.0) * temp2 * temp2) - ((1.0 / 24.0) * temp2)));   // (-1, -1, -1)
#else
  f_new[ f_id[i*15+7]+is_current*15 ] = f[7] += omega * (f_neq[7] = f[7] - ((temp1 + (1.0 / 16.0) * temp2 * temp2) + ((1.0 / 24.0) * temp2)));   // (+1, +1, +1)
  f_new[ f_id[i*15+8]+is_current*15 ] = f[8] += omega * (f_neq[8] = f[8] - ((temp1 + (1.0 / 16.0) * temp2 * temp2) - ((1.0 / 24.0) * temp2)));   // (-1, -1, -1)
#endif
  
  temp2 = (*v_x + *v_y) - *v_z;
  
#ifndef TD
  f_new[ f_id[i*15+ 9]               ] = f[ 9] += omega * (f_neq[ 9] = f[ 9] - ((temp1 + (1.0 / 16.0) * temp2 * temp2) + ((1.0 / 24.0) * temp2)));   // (+1, +1, -1)
  f_new[ f_id[i*15+10]               ] = f[10] += omega * (f_neq[10] = f[10] - ((temp1 + (1.0 / 16.0) * temp2 * temp2) - ((1.0 / 24.0) * temp2)));   // (-1, -1, +1)
#else
  f_new[ f_id[i*15+ 9]+is_current*15 ] = f[ 9] += omega * (f_neq[ 9] = f[ 9] - ((temp1 + (1.0 / 16.0) * temp2 * temp2) + ((1.0 / 24.0) * temp2)));   // (+1, +1, -1)
  f_new[ f_id[i*15+10]+is_current*15 ] = f[10] += omega * (f_neq[10] = f[10] - ((temp1 + (1.0 / 16.0) * temp2 * temp2) - ((1.0 / 24.0) * temp2)));   // (-1, -1, +1)
#endif
  
  temp2 = (*v_x - *v_y) + *v_z;
  
#ifndef TD
  f_new[ f_id[i*15+11]               ] = f[11] += omega * (f_neq[11] = f[11] - ((temp1 + (1.0 / 16.0) * temp2 * temp2) + ((1.0 / 24.0) * temp2)));   // (+1, -1, +1)
  f_new[ f_id[i*15+12]               ] = f[12] += omega * (f_neq[12] = f[12] - ((temp1 + (1.0 / 16.0) * temp2 * temp2) - ((1.0 / 24.0) * temp2)));   // (-1, +1, -1)
  #else
  f_new[ f_id[i*15+11]+is_current*15 ] = f[11] += omega * (f_neq[11] = f[11] - ((temp1 + (1.0 / 16.0) * temp2 * temp2) + ((1.0 / 24.0) * temp2)));   // (+1, -1, +1)
  f_new[ f_id[i*15+12]+is_current*15 ] = f[12] += omega * (f_neq[12] = f[12] - ((temp1 + (1.0 / 16.0) * temp2 * temp2) - ((1.0 / 24.0) * temp2)));   // (-1, +1, -1)
#endif
  
  temp2 = (*v_x - *v_y) - *v_z;	 
  
#ifndef TD
  f_new[ f_id[i*15+13]               ] = f[13] += omega * (f_neq[13] = f[13] - ((temp1 + (1.0 / 16.0) * temp2 * temp2) + ((1.0 / 24.0) * temp2)));   // (+1, -1, -1)
  f_new[ f_id[i*15+14]               ] = f[14] += omega * (f_neq[14] = f[14] - ((temp1 + (1.0 / 16.0) * temp2 * temp2) - ((1.0 / 24.0) * temp2)));   // (-1, +1, +1)
#else
  f_new[ f_id[i*15+13]+is_current*15 ] = f[13] += omega * (f_neq[13] = f[13] - ((temp1 + (1.0 / 16.0) * temp2 * temp2) + ((1.0 / 24.0) * temp2)));   // (+1, -1, -1)
  f_new[ f_id[i*15+14]+is_current*15 ] = f[14] += omega * (f_neq[14] = f[14] - ((temp1 + (1.0 / 16.0) * temp2 * temp2) - ((1.0 / 24.0) * temp2)));   // (-1, +1, +1)
#endif
  
#ifndef BENCH
  for (int l = 0; l < 15; l++)
    {
      if (f[l] < 0.) is_unstable = 1;
    }
#endif
}


void lbmInnerCollisionSIMD0 (double omega, int i,
			     double density[], double v_x[], double v_y[], double v_z[],
			     double f_neq[])
{
  double f_vec[15][SIMD_SIZE];
  double f_neq_vec[15][SIMD_SIZE];
  double v_xx, v_yy, v_zz;
  double temp1, temp2;
  
  int j, l;
  
  
  for (j = 0; j < SIMD_SIZE; j++)
    {
      for (l = 0; l < 15; l++)
	{
#ifndef TD
	  f_vec[l][j] = f_old[ (i+j)*15+l ];
#else
	  f_vec[l][j] = f_old[ (i+j)*30+is_current*15+l ];
#endif
	}
    }
#pragma ivdep
  for (j = 0; j < SIMD_SIZE; j++)
    {

      v_x[j] = f_vec[1][j] + (f_vec[7][j] + f_vec[9][j]) + (f_vec[11][j] + f_vec[13][j]);
      v_y[j] = f_vec[3][j] + (f_vec[12][j] + f_vec[14][j]);
      v_z[j] = f_vec[5][j] + f_vec[10][j];
      
      density[j] = f_vec[0][j] + (f_vec[2][j] + f_vec[4][j]) + (f_vec[6][j] + f_vec[8][j]) + v_x[j] + v_y[j] + v_z[j];
      
      v_x[j] -= f_vec[2][j] + (f_vec[8][j] + f_vec[10][j]) + (f_vec[12][j] + f_vec[14][j]);
      v_y[j] += (f_vec[7][j] + f_vec[9][j]) - (f_vec[4][j] + (f_vec[8][j] + f_vec[10][j]) + (f_vec[11][j] + f_vec[13][j]));
      v_z[j] += f_vec[7][j] + f_vec[11][j] + f_vec[14][j] - ((f_vec[6][j] + f_vec[8][j]) + f_vec[9][j] + f_vec[12][j] + f_vec[13][j]);
      
      v_xx = v_x[j] * v_x[j];
      v_yy = v_y[j] * v_y[j];
      v_zz = v_z[j] * v_z[j];
      
      temp1 = (1.0 / 8.0) * density[j];
      
      temp2 = v_xx + v_yy + v_zz;
      
      f_vec[0][j] += omega * (f_neq_vec[0][j] = f_vec[0][j] - (temp1 - (1.0 / 3.0) * temp2));
      
      temp1 -= (1.0 / 6.0) * temp2;
      
      f_vec[1][j] += omega * (f_neq_vec[1][j] = f_vec[1][j] - ((temp1 + 0.5 * v_xx) + ((1.0 / 3.0) * v_x[j])));   // (+1, 0, 0)
      f_vec[2][j] += omega * (f_neq_vec[2][j] = f_vec[2][j] - ((temp1 + 0.5 * v_xx) - ((1.0 / 3.0) * v_x[j])));   // (+1, 0, 0)
      
      f_vec[3][j] += omega * (f_neq_vec[3][j] = f_vec[3][j] - ((temp1 + 0.5 * v_yy) + ((1.0 / 3.0) * v_y[j])));   // (0, +1, 0)
      f_vec[4][j] += omega * (f_neq_vec[4][j] = f_vec[4][j] - ((temp1 + 0.5 * v_yy) - ((1.0 / 3.0) * v_y[j])));   // (0, -1, 0)
      
      f_vec[5][j] += omega * (f_neq_vec[5][j] = f_vec[5][j] - ((temp1 + 0.5 * v_zz) + ((1.0 / 3.0) * v_z[j])));   // (0, 0, +1)
      f_vec[6][j] += omega * (f_neq_vec[6][j] = f_vec[6][j] - ((temp1 + 0.5 * v_zz) - ((1.0 / 3.0) * v_z[j])));   // (0, 0, -1)
      
      temp1 *= (1.0 / 8.0);
      
      temp2 = (v_x[j] + v_y[j]) + v_z[j];
      
      f_vec[7][ j] += omega * (f_neq_vec[ 7][j] = f_vec[7][ j] - ((temp1 + (1.0 / 16.0) * temp2 * temp2) + ((1.0 / 24.0) * temp2)));   // (+1, +1, +1)
      f_vec[8][ j] += omega * (f_neq_vec[ 8][j] = f_vec[8][ j] - ((temp1 + (1.0 / 16.0) * temp2 * temp2) - ((1.0 / 24.0) * temp2)));   // (-1, -1, -1)
      
      temp2 = (v_x[j] + v_y[j]) - v_z[j];
      
      f_vec[9][ j] += omega * (f_neq_vec[ 9][j] = f_vec[9][ j] - ((temp1 + (1.0 / 16.0) * temp2 * temp2) + ((1.0 / 24.0) * temp2)));   // (+1, +1, -1)
      f_vec[10][j] += omega * (f_neq_vec[10][j] = f_vec[10][j] - ((temp1 + (1.0 / 16.0) * temp2 * temp2) - ((1.0 / 24.0) * temp2)));   // (-1, -1, +1)
      
      temp2 = (v_x[j] - v_y[j]) + v_z[j];
      
      f_vec[11][j] += omega * (f_neq_vec[11][j] = f_vec[11][j] - ((temp1 + (1.0 / 16.0) * temp2 * temp2) + ((1.0 / 24.0) * temp2)));   // (+1, -1, +1)
      f_vec[12][j] += omega * (f_neq_vec[12][j] = f_vec[12][j] - ((temp1 + (1.0 / 16.0) * temp2 * temp2) - ((1.0 / 24.0) * temp2)));   // (-1, +1, -1)
      
      temp2 = (v_x[j] - v_y[j]) - v_z[j];
      
      f_vec[13][j] += omega * (f_neq_vec[13][j] = f_vec[13][j] - ((temp1 + (1.0 / 16.0) * temp2 * temp2) + ((1.0 / 24.0) * temp2)));   // (+1, -1, -1)
      f_vec[14][j] += omega * (f_neq_vec[14][j] = f_vec[14][j] - ((temp1 + (1.0 / 16.0) * temp2 * temp2) - ((1.0 / 24.0) * temp2)));   // (-1, +1, +1)
    }
  for (j = 0; j < SIMD_SIZE; j++)
    {
      for (l = 0; l < 15; l++)
        {
#ifndef BENCH
	  if (f_vec[l][j] < 0.) is_unstable = 1;
#endif
	  
#ifndef TD
	  f_new[ f_id[(i+j)*15+l] ] = f_vec[l][j];
#else
	  f_new[ f_id[(i+j)*15+l]+is_current*15 ] = f_vec[l][j];
#endif 
          f_neq[ j*15+l ] = f_neq_vec[l][j];
        }
    }
}


void lbmInterCollisionSIMD0 (double omega, int i,
			     double density[], double v_x[], double v_y[], double v_z[],
			     double f_neq[])
{
  double f_vec[15][SIMD_SIZE];
  double f_neq_vec[15][SIMD_SIZE];
  double v_xx, v_yy, v_zz;
  double temp1, temp2;
  
  int j, l;
  
  
  for (j = 0; j < SIMD_SIZE; j++)
    {
      for (l = 0; l < 15; l++)
	{
#ifndef TD
	  f_vec[l][j] = f_old[ (i+j)*15+l ];
#else
	  f_vec[l][j] = f_old[ (i+j)*30+is_current*15+l ];
#endif
	}
    }
  for (j = 0; j < SIMD_SIZE; j++)
    {
      v_x[j] = f_vec[1][j] + (f_vec[7][j] + f_vec[9][j]) + (f_vec[11][j] + f_vec[13][j]);
      v_y[j] = f_vec[3][j] + (f_vec[12][j] + f_vec[14][j]);
      v_z[j] = f_vec[5][j] + f_vec[10][j];
      
      density[j] = f_vec[0][j] + (f_vec[2][j] + f_vec[4][j]) + (f_vec[6][j] + f_vec[8][j]) + v_x[j] + v_y[j] + v_z[j];
      
      v_x[j] -= f_vec[2][j] + (f_vec[8][j] + f_vec[10][j]) + (f_vec[12][j] + f_vec[14][j]);
      v_y[j] += (f_vec[7][j] + f_vec[9][j]) - (f_vec[4][j] + (f_vec[8][j] + f_vec[10][j]) + (f_vec[11][j] + f_vec[13][j]));
      v_z[j] += f_vec[7][j] + f_vec[11][j] + f_vec[14][j] - ((f_vec[6][j] + f_vec[8][j]) + f_vec[9][j] + f_vec[12][j] + f_vec[13][j]);
      
      v_xx = v_x[j] * v_x[j];
      v_yy = v_y[j] * v_y[j];
      v_zz = v_z[j] * v_z[j];
      
      temp1 = (1.0 / 8.0) * density[j];
      
      temp2 = v_xx + v_yy + v_zz;
      
      f_vec[0][j] += omega * (f_neq_vec[0][j] = f_vec[0][j] - (temp1 - (1.0 / 3.0) * temp2));
      
      temp1 -= (1.0 / 6.0) * temp2;
      
      f_vec[1][j] += omega * (f_neq_vec[1][j] = f_vec[1][j] - ((temp1 + 0.5 * v_xx) + ((1.0 / 3.0) * v_x[j])));   // (+1, 0, 0)
      f_vec[2][j] += omega * (f_neq_vec[2][j] = f_vec[2][j] - ((temp1 + 0.5 * v_xx) - ((1.0 / 3.0) * v_x[j])));   // (+1, 0, 0)
      
      f_vec[3][j] += omega * (f_neq_vec[3][j] = f_vec[3][j] - ((temp1 + 0.5 * v_yy) + ((1.0 / 3.0) * v_y[j])));   // (0, +1, 0)
      f_vec[4][j] += omega * (f_neq_vec[4][j] = f_vec[4][j] - ((temp1 + 0.5 * v_yy) - ((1.0 / 3.0) * v_y[j])));   // (0, -1, 0)
      
      f_vec[5][j] += omega * (f_neq_vec[5][j] = f_vec[5][j] - ((temp1 + 0.5 * v_zz) + ((1.0 / 3.0) * v_z[j])));   // (0, 0, +1)
      f_vec[6][j] += omega * (f_neq_vec[6][j] = f_vec[6][j] - ((temp1 + 0.5 * v_zz) - ((1.0 / 3.0) * v_z[j])));   // (0, 0, -1)
      
      temp1 *= (1.0 / 8.0);
      
      temp2 = (v_x[j] + v_y[j]) + v_z[j];
      
      f_vec[7][ j] += omega * (f_neq_vec[ 7][j] = f_vec[7][ j] - ((temp1 + (1.0 / 16.0) * temp2 * temp2) + ((1.0 / 24.0) * temp2)));   // (+1, +1, +1)
      f_vec[8][ j] += omega * (f_neq_vec[ 8][j] = f_vec[8][ j] - ((temp1 + (1.0 / 16.0) * temp2 * temp2) - ((1.0 / 24.0) * temp2)));   // (-1, -1, -1)
      
      temp2 = (v_x[j] + v_y[j]) - v_z[j];
      
      f_vec[9][ j] += omega * (f_neq_vec[ 9][j] = f_vec[9][ j] - ((temp1 + (1.0 / 16.0) * temp2 * temp2) + ((1.0 / 24.0) * temp2)));   // (+1, +1, -1)
      f_vec[10][j] += omega * (f_neq_vec[10][j] = f_vec[10][j] - ((temp1 + (1.0 / 16.0) * temp2 * temp2) - ((1.0 / 24.0) * temp2)));   // (-1, -1, +1)
      
      temp2 = (v_x[j] - v_y[j]) + v_z[j];
      
      f_vec[11][j] += omega * (f_neq_vec[11][j] = f_vec[11][j] - ((temp1 + (1.0 / 16.0) * temp2 * temp2) + ((1.0 / 24.0) * temp2)));   // (+1, -1, +1)
      f_vec[12][j] += omega * (f_neq_vec[12][j] = f_vec[12][j] - ((temp1 + (1.0 / 16.0) * temp2 * temp2) - ((1.0 / 24.0) * temp2)));   // (-1, +1, -1)
      
      temp2 = (v_x[j] - v_y[j]) - v_z[j];
      
      f_vec[13][j] += omega * (f_neq_vec[13][j] = f_vec[13][j] - ((temp1 + (1.0 / 16.0) * temp2 * temp2) + ((1.0 / 24.0) * temp2)));   // (+1, -1, -1)
      f_vec[14][j] += omega * (f_neq_vec[14][j] = f_vec[14][j] - ((temp1 + (1.0 / 16.0) * temp2 * temp2) - ((1.0 / 24.0) * temp2)));   // (-1, +1, +1)
    }
  for (j = 0; j < SIMD_SIZE; j++)
    {
      for (l = 0; l < 15; l++)
        {
#ifndef BENCH
	  if (f_vec[l][j] < 0.) is_unstable = 1;
#endif
	  
#ifndef TD
	  f_new[ f_id[(i+j)*15+l] ] = f_old[ (i+j)*15+l ] = f_vec[l][j];
#else
	  f_new[ f_id[(i+j)*15+l]+is_current*15 ] = f_old[ (i+j)*30+is_current*15+l ] = f_vec[l][j];
#endif 
          f_neq[ j*15+l ] = f_neq_vec[l][j];
        }
    }
}


void lbmCollision1 (double omega, int i,
		    double *density, double *v_x, double *v_y, double *v_z,
		    double f_neq[])
{
  double *f;
  double temp;
  
  int l;
  
  
#ifndef TD
  f = &f_old[ i*15 ];
#else
  f = &f_old[ i*30+is_current*15 ];
#endif
  
  for (l = 0; l < 15; l++)
    {
      f_neq[l] = f[l];
    }
  *v_x = *v_y = *v_z = 0.F;
  
  *density = 0.;
  
  for (l = 0; l < 15; l++) *density += f[l];
  
  temp = (1.0 / 8.0) * *density;

  for (l = 0; l < 7; l++) f[l] = temp;
  
  temp *= (1.0 / 8.0);
  
  for (l = 7; l < 15; l++) f[l] = temp;
  
  for (l = 0; l < 15; l++)
    {
#ifndef BENCH
      if (f[l] < 0.)  is_unstable = 1;
#endif
      
#ifndef TD
      f_neq[l] -= (f_new[ f_id[i*15+l] ] = f[l]);
#else
      f_neq[l] -= (f_new[ f_id[i*15+l]+is_current*15 ] = f[l]);
#endif
    }
}


void lbmCollisionSIMD1 (double omega, int i,
			double density[], double v_x[], double v_y[], double v_z[],
			double f_neq[])
{
  double *f;
  double temp;
  
  int j, l;
  
  
  for (j = 0; j < SIMD_SIZE; j++)
    {
#ifndef TD
      f = &f_old[ (i+j)*15 ];
#else
      f = &f_old[ (i+j)*30+is_current*15 ];
#endif
      for (l = 0; l < 15; l++)
	{
	  f_neq[ j*15+l ] = f[l];
	}
      v_x[j] = v_y[j] = v_z[j] = 0.F;
      
      density[j] = 0.;
      
      for (l = 0; l < 15; l++) density[j] += f[l];
      
      temp = (1.0 / 8.0) * density[j];
      
      for (l = 0; l < 7; l++) f[l] = temp;
      
      temp *= (1.0 / 8.0);
      
      for (l = 7; l < 15; l++) f[l] = temp;
      
      for (l = 0; l < 15; l++)
	{
#ifndef BENCH
	  if (f[l] < 0.) is_unstable = 1;
#endif
	  
#ifndef TD
	  f_neq[ j*15+l ] -= (f_new[ f_id[(i+j)*15+l] ] = f[l]);
#else
	  f_neq[ j*15+l ] -= (f_new[ f_id[(i+j)*15+l]+is_current*15 ] = f[l]);
#endif
	}
    }
}


void lbmCollision2 (double omega, int i,
		    double *density, double *v_x, double *v_y, double *v_z,
		    double f_neq[])
{
  double *f;
  double dummy_density;
  
  unsigned int boundary_id, l;
  
  
#ifndef TD
  f = &f_old[ i*15 ];
#else
  f = &f_old[ i*30+is_current*15 ];
#endif
  
  for (l = 0; l < 15; l++)
    {
      f_neq[l] = f[l];
    }
  boundary_id = (net_site_data[ i ] & BOUNDARY_ID_MASK) >> BOUNDARY_ID_SHIFT;
  
  *density = inlet_density[ boundary_id ];
  
  lbmDensityAndVelocity (f, &dummy_density, v_x, v_y, v_z);
  lbmFeq (*density, *v_x, *v_y, *v_z, f);
  
  for (l = 0; l < 15; l++)
    {
#ifndef BENCH
      if (f[l] < 0.) is_unstable = 1;
#endif
      
#ifndef TD
      f_new[ f_id[i*15+l] ] = f[l];
#else
      f_new[ f_id[i*15+l]+is_current*15 ] = f[l];
#endif
      f_neq[l] -= f[l];
    }
}


void lbmCollisionSIMD2 (double omega, int i,
			double density[], double v_x[], double v_y[], double v_z[],
			double f_neq[])
{
  double f_vec[15][SIMD_SIZE];
  double v_xx, v_yy, v_zz;
  double temp1, temp2;
  
  int j, l;
  
  unsigned int boundary_id;
  
  
  for (j = 0; j < SIMD_SIZE; j++)
    {
      for (l = 0; l < 15; l++)
	{
#ifndef TD
	  f_neq[ j*15+l ] = f_vec[l][j] = f_old[ (i+j)*15+l ];
#else
	  f_neq[ j*15+l ] = f_vec[l][j] = f_old[ (i+j)*30+is_current*15+l ];
#endif
	}
      boundary_id = (net_site_data[ i+j ] & BOUNDARY_ID_MASK) >> BOUNDARY_ID_SHIFT;
      density[j] = inlet_density[ boundary_id ];
    }
  for (j = 0; j < SIMD_SIZE; j++)
    {
      v_x[j] = f_vec[1][j] + (f_vec[7][j] + f_vec[9][j]) + (f_vec[11][j] + f_vec[13][j]);
      v_y[j] = f_vec[3][j] + (f_vec[12][j] + f_vec[14][j]);
      v_z[j] = f_vec[5][j] + f_vec[10][j];
      
      v_x[j] -= f_vec[2][j] + (f_vec[8][j] + f_vec[10][j]) + (f_vec[12][j] + f_vec[14][j]);
      v_y[j] += (f_vec[7][j] + f_vec[9][j]) - (f_vec[4][j] + (f_vec[8][j] + f_vec[10][j]) + (f_vec[11][j] + f_vec[13][j]));
      v_z[j] += f_vec[7][j] + f_vec[11][j] + f_vec[14][j] - ((f_vec[6][j] + f_vec[8][j]) + f_vec[9][j] + f_vec[12][j] + f_vec[13][j]);
      
      v_xx = v_x[j] * v_x[j];
      v_yy = v_y[j] * v_y[j];
      v_zz = v_z[j] * v_z[j];
      
      temp1 = (1.0 / 8.0) * density[j];
      
      temp2 = v_xx + v_yy + v_zz;
      
      f_vec[0][j] = temp1 - (1.0 / 3.0) * temp2;
      
      temp1 -= (1.0 / 6.0) * temp2;
      
      f_vec[1][j] = (temp1 + 0.5 * v_xx) + ((1.0 / 3.0) * v_x[j]);   // (+1, 0, 0)
      f_vec[2][j] = (temp1 + 0.5 * v_xx) - ((1.0 / 3.0) * v_x[j]);   // (+1, 0, 0)
					    
      f_vec[3][j] = (temp1 + 0.5 * v_yy) + ((1.0 / 3.0) * v_y[j]);   // (0, +1, 0)
      f_vec[4][j] = (temp1 + 0.5 * v_yy) - ((1.0 / 3.0) * v_y[j]);   // (0, -1, 0)
      
      f_vec[5][j] = (temp1 + 0.5 * v_zz) + ((1.0 / 3.0) * v_z[j]);   // (0, 0, +1)
      f_vec[6][j] = (temp1 + 0.5 * v_zz) - ((1.0 / 3.0) * v_z[j]);   // (0, 0, -1)
      
      temp1 *= (1.0 / 8.0);
      
      temp2 = (v_x[j] + v_y[j]) + v_z[j];
      
      f_vec[7][ j] = (temp1 + (1.0 / 16.0) * temp2 * temp2) + ((1.0 / 24.0) * temp2);   // (+1, +1, +1)
      f_vec[8][ j] = (temp1 + (1.0 / 16.0) * temp2 * temp2) - ((1.0 / 24.0) * temp2);   // (-1, -1, -1)
      
      temp2 = (v_x[j] + v_y[j]) - v_z[j];		     
      
      f_vec[9][ j] = (temp1 + (1.0 / 16.0) * temp2 * temp2) + ((1.0 / 24.0) * temp2);   // (+1, +1, -1)
      f_vec[10][j] = (temp1 + (1.0 / 16.0) * temp2 * temp2) - ((1.0 / 24.0) * temp2);   // (-1, -1, +1)
      
      temp2 = (v_x[j] - v_y[j]) + v_z[j];	     
      
      f_vec[11][j] = (temp1 + (1.0 / 16.0) * temp2 * temp2) + ((1.0 / 24.0) * temp2);   // (+1, -1, +1)
      f_vec[12][j] = (temp1 + (1.0 / 16.0) * temp2 * temp2) - ((1.0 / 24.0) * temp2);   // (-1, +1, -1)
      
      temp2 = (v_x[j] - v_y[j]) - v_z[j];     
      
      f_vec[13][j] = (temp1 + (1.0 / 16.0) * temp2 * temp2) + ((1.0 / 24.0) * temp2);   // (+1, -1, -1)
      f_vec[14][j] = (temp1 + (1.0 / 16.0) * temp2 * temp2) - ((1.0 / 24.0) * temp2);   // (-1, +1, +1)
    }
  for (j = 0; j < SIMD_SIZE; j++)
    {
      for (l = 0; l < 15; l++)
        {
#ifndef BENCH
	  if (f_vec[l][j] < 0.) is_unstable = 1;
#endif
	  
#ifndef TD
          f_new[ f_id[(i+j)*15+l] ] = f_old[ (i+j)*15+l ] = f_vec[l][j];
#else
	  f_new[ f_id[(i+j)*15+l]+is_current*15 ] = f_old[ (i+j)*30+is_current*15+l ] = f_vec[l][j];
#endif
	  f_neq[ j*15+l ] -= f_vec[l][j];
        }
    }
}


void lbmCollision3 (double omega, int i,
		    double *density, double *v_x, double *v_y, double *v_z,
		    double f_neq[])
{
  double *f;
  double dummy_density;
  
  unsigned int boundary_id, l;
  
  
#ifndef TD
  f = &f_old[ i*15 ];
#else
  f = &f_old[ i*30+is_current*15 ];
#endif
  
  for (l = 0; l < 15; l++)
    {
      f_neq[l] = f[l];
    }
  boundary_id = (net_site_data[ i ] & BOUNDARY_ID_MASK) >> BOUNDARY_ID_SHIFT;
  
  *density = outlet_density[ boundary_id ];
  
  lbmDensityAndVelocity (f, &dummy_density, v_x, v_y, v_z);
  lbmFeq (*density, *v_x, *v_y, *v_z, f);
  
  for (l = 0; l < 15; l++)
    {
#ifndef BENCH
      if (f[l] < 0.) is_unstable = 1;
#endif
      
#ifndef TD
      f_new[ f_id[i*15+l] ] = f[l];
#else
      f_new[ f_id[i*15+l]+is_current*15 ] = f[l];
#endif
      f_neq[l] -= f[l];
    }
}


void lbmCollisionSIMD3 (double omega, int i,
			double density[], double v_x[], double v_y[], double v_z[],
			double f_neq[])
{
  double f_vec[15][SIMD_SIZE];
  double v_xx, v_yy, v_zz;
  double temp1, temp2;
  
  int j, l;
  
  unsigned int boundary_id;
  
  
  for (j = 0; j < SIMD_SIZE; j++)
    {
      for (l = 0; l < 15; l++)
	{
#ifndef TD
	  f_neq[ j*15+l ] = f_vec[l][j] = f_old[ (i+j)*15+l ];
#else
	  f_neq[ j*15+l ] = f_vec[l][j] = f_old[ (i+j)*30+is_current*15+l ];
#endif
	}
      boundary_id = (net_site_data[ i+j ] & BOUNDARY_ID_MASK) >> BOUNDARY_ID_SHIFT;
      density[j] = outlet_density[ boundary_id ];
    }
  for (j = 0; j < SIMD_SIZE; j++)
    {
      v_x[j] = f_vec[1][j] + (f_vec[7][j] + f_vec[9][j]) + (f_vec[11][j] + f_vec[13][j]);
      v_y[j] = f_vec[3][j] + (f_vec[12][j] + f_vec[14][j]);
      v_z[j] = f_vec[5][j] + f_vec[10][j];
      
      v_x[j] -= f_vec[2][j] + (f_vec[8][j] + f_vec[10][j]) + (f_vec[12][j] + f_vec[14][j]);
      v_y[j] += (f_vec[7][j] + f_vec[9][j]) - (f_vec[4][j] + (f_vec[8][j] + f_vec[10][j]) + (f_vec[11][j] + f_vec[13][j]));
      v_z[j] += f_vec[7][j] + f_vec[11][j] + f_vec[14][j] - ((f_vec[6][j] + f_vec[8][j]) + f_vec[9][j] + f_vec[12][j] + f_vec[13][j]);
      
      v_xx = v_x[j] * v_x[j];
      v_yy = v_y[j] * v_y[j];
      v_zz = v_z[j] * v_z[j];
      
      temp1 = (1.0 / 8.0) * density[j];
      
      temp2 = v_xx + v_yy + v_zz;
      
      f_vec[0][j] = temp1 - (1.0 / 3.0) * temp2;
      
      temp1 -= (1.0 / 6.0) * temp2;
      
      f_vec[1][j] = (temp1 + 0.5 * v_xx) + ((1.0 / 3.0) * v_x[j]);   // (+1, 0, 0)
      f_vec[2][j] = (temp1 + 0.5 * v_xx) - ((1.0 / 3.0) * v_x[j]);   // (+1, 0, 0)
      
      f_vec[3][j] = (temp1 + 0.5 * v_yy) + ((1.0 / 3.0) * v_y[j]);   // (0, +1, 0)
      f_vec[4][j] = (temp1 + 0.5 * v_yy) - ((1.0 / 3.0) * v_y[j]);   // (0, -1, 0)
      
      f_vec[5][j] = (temp1 + 0.5 * v_zz) + ((1.0 / 3.0) * v_z[j]);   // (0, 0, +1)
      f_vec[6][j] = (temp1 + 0.5 * v_zz) - ((1.0 / 3.0) * v_z[j]);   // (0, 0, -1)
      
      temp1 *= (1.0 / 8.0);
      
      temp2 = (v_x[j] + v_y[j]) + v_z[j];
      
      f_vec[7][ j] = (temp1 + (1.0 / 16.0) * temp2 * temp2) + ((1.0 / 24.0) * temp2);   // (+1, +1, +1)
      f_vec[8][ j] = (temp1 + (1.0 / 16.0) * temp2 * temp2) - ((1.0 / 24.0) * temp2);   // (-1, -1, -1)
      
      temp2 = (v_x[j] + v_y[j]) - v_z[j];		     
      
      f_vec[9][ j] = (temp1 + (1.0 / 16.0) * temp2 * temp2) + ((1.0 / 24.0) * temp2);   // (+1, +1, -1)
      f_vec[10][j] = (temp1 + (1.0 / 16.0) * temp2 * temp2) - ((1.0 / 24.0) * temp2);   // (-1, -1, +1)
      
      temp2 = (v_x[j] - v_y[j]) + v_z[j];	     
      
      f_vec[11][j] = (temp1 + (1.0 / 16.0) * temp2 * temp2) + ((1.0 / 24.0) * temp2);   // (+1, -1, +1)
      f_vec[12][j] = (temp1 + (1.0 / 16.0) * temp2 * temp2) - ((1.0 / 24.0) * temp2);   // (-1, +1, -1)
      
      temp2 = (v_x[j] - v_y[j]) - v_z[j];     
      
      f_vec[13][j] = (temp1 + (1.0 / 16.0) * temp2 * temp2) + ((1.0 / 24.0) * temp2);   // (+1, -1, -1)
      f_vec[14][j] = (temp1 + (1.0 / 16.0) * temp2 * temp2) - ((1.0 / 24.0) * temp2);   // (-1, +1, +1)
    }
  for (j = 0; j < SIMD_SIZE; j++)
    {
      for (l = 0; l < 15; l++)
        {
#ifndef BENCH
	  if (f_vec[l][j] < 0.) is_unstable = 1;
#endif
	  
#ifndef TD
          f_new[ f_id[(i+j)*15+l] ] = f_old[ (i+j)*15+l ] = f_vec[l][j];
#else
	  f_new[ f_id[(i+j)*15+l]+is_current*15 ] = f_old[ (i+j)*30+is_current*15+l ] = f_vec[l][j];
#endif
	  f_neq[ j*15+l ] -= f_vec[l][j];
        }
    }
}


void lbmCollision4 (double omega, int i,
		    double *density, double *v_x, double *v_y, double *v_z,
		    double f_neq[])
{
  double *f;
  double temp;
  
  int l;
  
  unsigned int boundary_id;
  
  
#ifndef TD
  f = &f_old[ i*15 ];
#else
  f = &f_old[ i*30+is_current*15 ];
#endif
  
  for (l = 0; l < 15; l++)
    {
      f_neq[l] = f[l];
    }
  boundary_id = (net_site_data[ i ] & BOUNDARY_ID_MASK) >> BOUNDARY_ID_SHIFT;
  
  *density = inlet_density[ boundary_id ];
  
  *v_x = *v_y = *v_z = 0.F;
  
  temp = (1.0 / 8.0) * *density;
  
  for (l = 0; l < 7; l++) f[l] = temp;
  
  temp *= (1.0 / 8.0);
  
  for (l = 7; l < 15; l++) f[l] = temp;
  
  for (l = 0; l < 15; l++)
    {
#ifndef BENCH
      if (f[l] < 0.) is_unstable = 1;
#endif
      
#ifndef TD
      f_neq[l] -= (f_new[ f_id[i*15+l] ] = f[l]);
#else
      f_neq[l] -= (f_new[ f_id[i*15+l]+is_current*15 ] = f[l]);
#endif
    }
}


void lbmCollisionSIMD4 (double omega, int i,
			double density[], double v_x[], double v_y[], double v_z[],
			double f_neq[])
{
  double *f;
  double temp;
  
  int j, l;
  
  unsigned int boundary_id;
  
  
  for (j = 0; j < SIMD_SIZE; j++)
    {
#ifndef TD
      f = &f_old[ (i+j)*15 ];
#else
      f = &f_old[ (i+j)*30+is_current*15 ];
#endif
      for (l = 0; l < 15; l++)
	{
	  f_neq[ j*15+l ] = f[l];
	}
      boundary_id = (net_site_data[ i+j ] & BOUNDARY_ID_MASK) >> BOUNDARY_ID_SHIFT;
      
      v_x[j] = v_y[j] = v_z[j] = 0.F;
      
      density[j] = inlet_density[ boundary_id ];
      
      temp = (1.0 / 8.0) * density[j];
      
      for (l = 0; l < 7; l++) f[l] = temp;
      
      temp *= (1.0 / 8.0);
      
      for (l = 7; l < 15; l++) f[l] = temp;
      
      for (l = 0; l < 15; l++)
	{
#ifndef BENCH
	  if (f[l] < 0.) is_unstable = 1;
#endif
	  
#ifndef TD
	  f_neq[ j*15+l ] -= (f_new[ f_id[(i+j)*15+l] ] = f[l]);
#else
	  f_neq[ j*15+l ] -= (f_new[ f_id[(i+j)*15+l]+is_current*15 ] = f[l]);
#endif
	}
    }
}


void lbmCollision5 (double omega, int i,
		    double *density, double *v_x, double *v_y, double *v_z,
		    double f_neq[])
{
  double *f;
  double temp;
  
  int l;
  
  unsigned int boundary_id;
  
  
#ifndef TD
  f = &f_old[ i*15 ];
#else
  f = &f_old[ i*30+is_current*15 ];
#endif
  
  for (l = 0; l < 15; l++)
    {
      f_neq[l] = f[l];
    }
  boundary_id = (net_site_data[ i ] & BOUNDARY_ID_MASK) >> BOUNDARY_ID_SHIFT;
  
  *density = outlet_density[ boundary_id ];
  
  *v_x = *v_y = *v_z = 0.F;
  
  temp = (1.0 / 8.0) * *density;
  
  for (l = 0; l < 7; l++) f[l] = temp;
  
  temp *= (1.0 / 8.0);
  
  for (l = 7; l < 15; l++) f[l] = temp;
  
  for (l = 0; l < 15; l++)
    {
#ifndef BENCH
      if (f[l] < 0.) is_unstable = 1;
#endif
      
#ifndef TD
      f_neq[l] -= (f_new[ f_id[i*15+l] ] = f[l]);
#else
      f_neq[l] -= (f_new[ f_id[i*15+l]+is_current*15 ] = f[l]);
#endif
    }
}


void lbmCollisionSIMD5 (double omega, int i,
			double density[], double v_x[], double v_y[], double v_z[],
			double f_neq[])
{
  double *f;
  double temp;
  
  int j, l;
  
  unsigned int boundary_id;
  
  
  for (j = 0; j < SIMD_SIZE; j++)
    {
#ifndef TD
      f = &f_old[ (i+j)*15 ];
#else
      f = &f_old[ (i+j)*30+is_current*15 ];
#endif
      for (l = 0; l < 15; l++)
	{
	  f_neq[ j*15+l ] = f[ l ];
	}
      boundary_id = (net_site_data[ i+j ] & BOUNDARY_ID_MASK) >> BOUNDARY_ID_SHIFT;
      
      v_x[j] = v_y[j] = v_z[j] = 0.F;
      
      density[j] = outlet_density[ boundary_id ];
      
      temp = (1.0 / 8.0) * density[j];
      
      for (l = 0; l < 7; l++) f[l] = temp;
      
      temp *= (1.0 / 8.0);
      
      for (l = 7; l < 15; l++) f[l] = temp;
      
      for (l = 0; l < 15; l++)
	{
#ifndef BENCH
	  if (f[l] < 0.) is_unstable = 1;
#endif
	  
#ifndef TD
	  f_neq[ j*15+l ] -= (f_new[ f_id[(i+j)*15+l] ] = f[l]);
#else
	  f_neq[ j*15+l ] -= (f_new[ f_id[(i+j)*15+l]+is_current*15 ] = f[l]);
#endif
	}
    }
}


void lbmVelocity (double f[], double *v_x, double *v_y, double *v_z)
{
  *v_x = f[1] + (f[ 7] + f[  9 ]) + (f[ 11 ] + f[ 13 ]);
  *v_y = f[3] + (f[ 12 ] + f[ 14 ]);
  *v_z = f[5] + f[ 10 ];
  
  *v_x -= (f[2] + f[ 8 ] + f[ 10 ] + (f[ 12 ] + f[ 14 ]));
  *v_y += (f[7] + f[ 9 ]) - ((f[4] + f[ 8 ] + f[ 10 ] + (f[ 11 ] + f[ 13 ])));
  *v_z += f[7] + f[ 11 ] + f[ 14 ] - (((f[6] + f[ 8 ]) + f[ 9 ] + f[ 12 ] + f[ 13 ]));
}


void lbmDensityAndVelocity (double f[], double *density, double *v_x, double *v_y, double *v_z)
{
  *v_x = f[1] + (f[7] + f[9]) + (f[11] + f[13]);
  *v_y = f[3] + (f[12] + f[14]);
  *v_z = f[5] + f[10];
  
  *density = f[0] + (f[2] + f[4]) + (f[6] + f[8]) + *v_x + *v_y + *v_z;
  
  *v_x -= (f[2] + f[8] + f[10] + (f[12] + f[14]));
  *v_y += (f[7] + f[9]) - ((f[4] + f[8] + f[10] + (f[11] + f[13])));
  *v_z += f[7] + f[11] + f[14] - (((f[6] + f[8]) + f[9] + f[12] + f[13]));
}


void lbmStress (double f[], double *stress)
{
  double sigma_xx_yy, sigma_yy_zz, sigma_xx_zz;
  double sigma_xy, sigma_xz, sigma_yz;
  double a, b;
  
  
  sigma_xx_yy = (f[1] + f[2]) - (f[3] + f[4]);
  sigma_yy_zz = (f[3] + f[4]) - (f[5] + f[6]);
  sigma_xx_zz = (f[1] + f[2]) - (f[5] + f[6]);
  
  sigma_xy = (f[7] + f[8]) + (f[9] + f[10]) - (f[11] + f[12]) - (f[13] + f[14]);
  sigma_xz = (f[7] + f[8]) - (f[9] + f[10]) + (f[11] + f[12]) - (f[13] + f[14]);
  sigma_yz = (f[7] + f[8]) - (f[9] + f[10]) - (f[11] + f[12]) + (f[13] + f[14]);
  
  a = sigma_xx_yy * sigma_xx_yy + sigma_yy_zz * sigma_yy_zz + sigma_xx_zz * sigma_xx_zz;
  b = sigma_xy * sigma_xy + sigma_xz * sigma_xz + sigma_yz * sigma_yz;
  
  *stress = lbm_stress_par * sqrt(a + 6.0 * b);
}


void lbmUpdateFlowField (int i, double density, double vx, double vy, double vz, double f_neq[])
{
  double velocity, stress;
  
  
  velocity = sqrt(vx * vx + vy * vy + vz * vz);
  
  lbmStress (f_neq, &stress);
  
  lbm_density_min = (density < lbm_density_min) ? density : lbm_density_min;
  lbm_density_max = (density > lbm_density_max) ? density : lbm_density_max;
  
  lbm_velocity_min = (velocity < lbm_velocity_min) ? velocity : lbm_velocity_min;
  lbm_velocity_max = (velocity > lbm_velocity_max) ? velocity : lbm_velocity_max;
  
  lbm_stress_min = (stress < lbm_stress_min) ? stress : lbm_stress_min;
  lbm_stress_max = (stress > lbm_stress_max) ? stress : lbm_stress_max;
  
  if (vis_flow_field_type == DENSITY)
    {
      *cluster_voxel[ i ] = (float)density;
    }
  else if (vis_flow_field_type == VELOCITY)
    {
      *cluster_voxel[ i ] = (float)velocity;
    }
  else
    {
      lbmStress (f_neq, &stress);
      
      *cluster_voxel[ i ] = (float)stress;
    }
}


int lbmCollisionType (unsigned int site_data)
{
  unsigned int boundary_type, boundary_config;
  
  int unknowns, i;
  
  
  if (site_data == FLUID_TYPE)
    {
      return 0;
    }
  boundary_type = site_data & SITE_TYPE_MASK;
  
  if (boundary_type == FLUID_TYPE)
    {
      return 1;
    }
  boundary_config = (site_data & BOUNDARY_CONFIG_MASK) >> BOUNDARY_CONFIG_SHIFT;
  
  unknowns = 0;
  
  for (i = 0; i < 14; i++)
    {
      if (!(boundary_config & (1U << i))) ++unknowns;
    }
  if (boundary_type == INLET_TYPE)
    {
      if (unknowns <= 5)
	{
	  return 2;
	}
      else
	{
	  return 4;
	}
    }
  else
    {
      if (unknowns <= 5)
	{
	  return 3;
	}
      else
	{
	  return 5;
	}
    }
}


void lbmCalculateBC (double f[], unsigned int site_data, double *density,
		     double *vx, double *vy, double *vz, double f_neq[])
{
  double dummy_density;
  double temp;
  
  int unknowns, i;
  
  unsigned int boundary_type, boundary_config, boundary_id;
  
  

  
  for (i = 0; i < 15; i++)
    {
      f_neq[ i ] = f[ i ];
    }
  boundary_type = site_data & SITE_TYPE_MASK;
  
  if (boundary_type == FLUID_TYPE)
    {
      *density = 0.;

      for (i = 0; i < 15; i++) *density += f[ i ];
      
      temp = *density * (1.0 / 8.0);
      
      for (i = 0; i < 7; i++) f[ i ] = temp;
      
      temp *= (1.0 / 8.0);
      
      for (i = 7; i < 15; i++) f[ i ] = temp;
      
      *vx = *vy = *vz = 0.F;
    }
  else
    {
      boundary_config = (site_data & BOUNDARY_CONFIG_MASK) >> BOUNDARY_CONFIG_SHIFT;
      boundary_id     = (site_data & BOUNDARY_ID_MASK)     >> BOUNDARY_ID_SHIFT;
      
      unknowns = 0;
      
      for (i = 0; i < 14; i++)
	{
	  if (!(boundary_config & (1U << i))) ++unknowns;
	}
      if (boundary_type == INLET_TYPE)
	{
	  *density = inlet_density[ boundary_id ];
	}
      else
	{
	  *density = outlet_density[ boundary_id ];
	}
      if (unknowns <= 5)
	{
	  lbmDensityAndVelocity (f, &dummy_density, vx, vy, vz);
	  lbmFeq (*density, *vx, *vy, *vz, f);
	}
      else
	{
	  temp = *density * (1.0 / 8.0);
	  
	  for (i = 0; i < 7; i++) f[ i ] = temp;
	  
	  temp *= (1.0 / 8.0);
	  
	  for (i = 7; i < 15; i++) f[ i ] = temp;
	  
	  *vx = *vy = *vz = 0.F;
	}
    }
  for (i = 0; i < 15; i++)
    {
      f_neq[ i ] -= f[ i ];
    }
}


void lbmInit (char *system_file_name, char *checkpoint_file_name,
	      LBM *lbm, Net *net)
{
  lbm->system_file_name     = system_file_name;
  lbm->checkpoint_file_name = checkpoint_file_name;
  
  lbmReadConfig (lbm, net);
  
  lbmInnerCollision[0] = lbmInnerCollision0;
  lbmInnerCollision[1] = lbmCollision1;
  lbmInnerCollision[2] = lbmCollision2;
  lbmInnerCollision[3] = lbmCollision3;
  lbmInnerCollision[4] = lbmCollision4;
  lbmInnerCollision[5] = lbmCollision5;
  
  lbmInterCollision[0] = lbmInterCollision0;
  lbmInterCollision[1] = lbmCollision1;
  lbmInterCollision[2] = lbmCollision2;
  lbmInterCollision[3] = lbmCollision3;
  lbmInterCollision[4] = lbmCollision4;
  lbmInterCollision[5] = lbmCollision5;
  
  lbmInnerCollisionSIMD[0] = lbmInnerCollisionSIMD0;
  lbmInnerCollisionSIMD[1] = lbmCollisionSIMD1;
  lbmInnerCollisionSIMD[2] = lbmCollisionSIMD2;
  lbmInnerCollisionSIMD[3] = lbmCollisionSIMD3;
  lbmInnerCollisionSIMD[4] = lbmCollisionSIMD4;
  lbmInnerCollisionSIMD[5] = lbmCollisionSIMD5;
  
  lbmInterCollisionSIMD[0] = lbmInterCollisionSIMD0;
  lbmInterCollisionSIMD[1] = lbmCollisionSIMD1;
  lbmInterCollisionSIMD[2] = lbmCollisionSIMD2;
  lbmInterCollisionSIMD[3] = lbmCollisionSIMD3;
  lbmInterCollisionSIMD[4] = lbmCollisionSIMD4;
  lbmInterCollisionSIMD[5] = lbmCollisionSIMD5;
}


void lbmSetInitialConditions (Net *net)
{
  double *f_old_p, *f_new_p, f_eq[15];
  
  int i, l;
  
  
#ifdef BENCH
  
  for (i = 0; i < net->my_sites; i++)
    {
      for (l = 0; l < 7; l++) f_eq[ l ] = (1.0 / 8.0);
      
      for (l = 7; l < 15; l++) f_eq[ l ] = (1.0 / 64.0);
#ifndef TD
      f_old_p = &f_old[ i*15 ];
      f_new_p = &f_new[ i*15 ];
      
      for (l = 0; l < 15; l++)
	{
	  f_new_p[ l ] = f_old_p[ l ] = f_eq[ l ];
	}
#else
      for (int cycle_id = 0; cycle_id < 2; cycle_id++)
	{
	  f_old_p = &f_old[ (i*2+cycle_id)*15 ];
	  f_new_p = &f_new[ (i*2+cycle_id)*15 ];
	  
	  for (l = 0; l < 15; l++)
	    {
	      f_new_p[ l ] = f_old_p[ l ] = f_eq[ l ];
	    }
	}
#endif
    }
  return;
#else   // BENCH
  
  double *d_p;
  double **nd_p_p;
  double density;
  double error, error_tot;
  double temp;
  
  int iters;
  int neighs;
  int m, n;
  
  unsigned int site_data, boundary_id;
  
  NeighProc *neigh_proc_p;
  
  
  for (i = 0; i < net->my_sites; i++)
    {
      site_data = net_site_data[ i ];
      
      if ((site_data & SITE_TYPE_MASK) == INLET_TYPE)
	{
	  boundary_id = (site_data & BOUNDARY_ID_MASK) >> BOUNDARY_ID_SHIFT;
	  
	  d[ i ] = inlet_density[ boundary_id ];
	}
      else if ((site_data & SITE_TYPE_MASK) == OUTLET_TYPE)
	{
	  boundary_id = (site_data & BOUNDARY_ID_MASK) >> BOUNDARY_ID_SHIFT;
	  
	  d[ i ] = outlet_density[ boundary_id ];
	}
    }
  d[ net->my_sites ] = -1.;
  
  error_tot = 1.e+30;
  iters = 0;
  
  while (error_tot > 1.e-3)
    {
      ++iters;
      
      error = 0.;
      
      for (m = 0; m < net->neigh_procs; m++)
	{
	  neigh_proc_p = &net->neigh_proc[ m ];
	  
	  for (n = 0; n < neigh_proc_p->fs; n++)
	    {
	      neigh_proc_p->f_to_send[ n ] = *neigh_proc_p->d_to_send_p[ n ];
	    }
#ifndef NOMPI
	  net->err = MPI_Isend (&neigh_proc_p->f_to_send[ 0 ],
				neigh_proc_p->fs, MPI_DOUBLE,
				neigh_proc_p->id, 10, MPI_COMM_WORLD,
				&net->req[ 0 ][ m ]);
	  
	  net->err = MPI_Irecv (&neigh_proc_p->f_to_recv[ 0 ],
				neigh_proc_p->fs, MPI_DOUBLE,
				neigh_proc_p->id, 10, MPI_COMM_WORLD,
				&net->req[ 0 ][ net->neigh_procs + m ]);
#endif
	}
      for (i = 0; i < net->my_inner_sites; i++)
	{
	  site_data = net_site_data[ i ];
	  
	  if ((site_data & SITE_TYPE_MASK) == INLET_TYPE ||
	      (site_data & SITE_TYPE_MASK) == OUTLET_TYPE)
	    {
	      continue;
	    }
	  d_p = &d[ i ];
	  
	  temp = *d_p;
	  *d_p = 0.;
	  neighs = 0;
	  
	  nd_p_p = &nd_p[ i*14 ];
	  
	  for (l = 0; l < 14; l++)
	    {
	      if (*nd_p_p[ l ] < 0.) continue;
	      
	      ++neighs;
	      *d_p += *nd_p_p[ l ];
	    }
	  *d_p /= (double)neighs;
	  
	  error = fmax(error, fabs(*d_p - temp) / fmax(1.e-30, *d_p));
	}
      for (m = 0; m < net->neigh_procs; m++)
      	{
#ifndef NOMPI
      	  net->err = MPI_Wait (&net->req[ 0 ][ m ], net->status);
      	  net->err = MPI_Wait (&net->req[ 0 ][ net->neigh_procs + m ], net->status);
#endif
      	}
      for (i = net->my_inner_sites; i < net->my_sites; i++)
	{
	  site_data = net_site_data[ i ];
	  
	  if ((site_data & SITE_TYPE_MASK) == INLET_TYPE ||
	      (site_data & SITE_TYPE_MASK) == OUTLET_TYPE)
	    {
	      continue;
	    }
	  d_p = &d[ i ];
	  
	  temp = *d_p;
	  *d_p = 0.;
	  neighs = 0;
	  
	  nd_p_p = &nd_p[ i*14 ];
	  
	  for (l = 0; l < 14; l++)
	    {
	      if (*nd_p_p[ l ] < 0.) continue;
	      
	      ++neighs;
	      *d_p += *nd_p_p[ l ];
	    }
	  *d_p /= (double)neighs;
	  
	  error = fmax(error, fabs(*d_p - temp) / fmax(1.e-30, *d_p));
	}
#ifndef NOMPI
      net->err = MPI_Allreduce (&error, &error_tot, 1,
				MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
      error_tot = error;
#endif
    }
  for (i = 0; i < net->my_sites; i++)
    {
      density = d[ i ];
      
      density *= (1.0 / 8.0);
      
      for (l = 0; l < 7; l++) f_eq[ l ] = density;
      
      density *= (1.0 / 8.0);
      
      for (l = 7; l < 15; l++) f_eq[ l ] = density;
#ifndef TD
      f_old_p = &f_old[ i*15 ];
      f_new_p = &f_new[ i*15 ];
      
      for (l = 0; l < 15; l++)
	{
	  f_new_p[ l ] = f_old_p[ l ] = f_eq[ l ];
	}
#else
      for (int cycle_id = 0; cycle_id < 2; cycle_id++)
	{
	  f_old_p = &f_old[ (i*2+cycle_id)*15 ];
	  f_new_p = &f_new[ (i*2+cycle_id)*15 ];
	  
	  for (l = 0; l < 15; l++)
	    {
	      f_new_p[ l ] = f_old_p[ l ] = f_eq[ l ];
	    }
	}
#endif
    }
  
  for (n = 0; n < net->neigh_procs; n++)
    {
      free(net->neigh_proc[ n ].d_to_send_p);
    }
  free (nd_p);
  nd_p = NULL;
  
  free (d);
  d = NULL;
#endif
}

/*
// old version

int lbmCycle (int write_checkpoint, int check_conv, int *is_converged, LBM *lbm, Net *net)
{
  // the entire simulation time step takes place through this function
  
  double f_eq[15], f_neq[15];
  double omega;
  double density;
  double vx, vy, vz;
  double *f_old_p;
  
#ifndef BENCH
  double sum1, sum2;
  double local_data[6];
  double global_data[6];
#endif // BENCH
  
  int unit_level;
  int i, l, m, n;
  int is_unstable;
  
  NeighProc *neigh_proc_p;
  
  
  is_unstable = 0;
  
  *is_converged = 0;
  
#ifndef BENCH
  sum1 = 0.0;
  sum2 = 0.0;
#endif // BENCH
  
  omega = lbm->omega;
  
  for (i = net->my_inner_sites; i < net->my_inner_sites + net->my_inter_sites; i++)
    {
      if (net_site_data[ i ] == FLUID_TYPE)
	{
	  lbmFeq (&f_old[ i*15 ], &density, &vx, &vy, &vz, f_eq);
	  
	  for (l = 0; l < 15; l++)
	    {
	      f_old[ i*15+l ] += omega * (f_neq[l] = f_old[ i*15+l ] - f_eq[l]);
	    }
	}
      else
	{
	  lbmCalculateBC (&f_old[ i*15 ], net_site_data[ i ], &density, &vx, &vy, &vz, f_neq);
	}
#ifndef BENCH
      sum1 += fabs(vel[ 3*i+0 ] - vx) + fabs(vel[ 3*i+1 ] - vy) + fabs(vel[ 3*i+2 ] - vz);
      sum2 += fabs(vx) + fabs(vy) + fabs(vz);
      
      vel[ 3*i+0 ] = vx;
      vel[ 3*i+1 ] = vy;
      vel[ 3*i+2 ] = vz;
#endif // BENCH
      
      lbmUpdateFlowField (i, density, vx, vy, vz, f_neq);
    }
  
  for (unit_level = 1; unit_level >= 0; unit_level--)
    {
      for (m = 0; m < net->neigh_procs; m++)
	{
	  neigh_proc_p = &net->neigh_proc[ m ];
	  
	  if ((unit_level == 1 && net->machine_id[ neigh_proc_p->id ] == net->machine_id[ net->id ]) ||
	      (unit_level == 0 && net->machine_id[ neigh_proc_p->id ] != net->machine_id[ net->id ]))
	    {
	      continue;
	    }
	  for (n = 0; n < neigh_proc_p->fs; n++)
	    {
	      neigh_proc_p->f_to_send[ n ] = f_old[ neigh_proc_p->f_send_id[n] ];
	    }
#ifndef NOMPI
	  net->err = MPI_Isend (&neigh_proc_p->f_to_send[ 0 ],
				neigh_proc_p->fs, MPI_DOUBLE,
				neigh_proc_p->id, 10, MPI_COMM_WORLD,
				&net->req[ 0 ][ m ]);
	  
	  net->err = MPI_Irecv (&neigh_proc_p->f_to_recv[ 0 ],
				neigh_proc_p->fs, MPI_DOUBLE,
				neigh_proc_p->id, 10, MPI_COMM_WORLD,
				&net->req[ 0 ][ net->neigh_procs + m ]);
#endif
	}
    }
  
  for (i = 0; i < net->my_inner_sites; i++)
    {
      if (net_site_data[ i ] == FLUID_TYPE)
	{
	  lbmFeq (&f_old[ i*15 ], &density, &vx, &vy, &vz, f_eq);
	  
	  for (l = 0; l < 15; l++)
	    {
	      f_old[ i*15+l ] += omega * (f_neq[l] = f_old[ i*15+l ] - f_eq[l]);
	    }
	}
      else
	{
	  lbmCalculateBC (&f_old[ i*15 ], net_site_data[ i ], &density, &vx, &vy, &vz, f_neq);
	}
      for (l = 0; l < 15; l++)
	{
#ifndef BENCH
	  if (f_old[ i*15+l ] < 0.) is_unstable = 1;
#endif // BENCH
	  
	  f_new[ f_id[i*15+l] ] = f_old[ i*15+l ];
	}
#ifndef BENCH
      sum1 += fabs(vel[ 3*i+0 ] - vx) + fabs(vel[ 3*i+1 ] - vy) + fabs(vel[ 3*i+2 ] - vz);
      sum2 += fabs(vx) + fabs(vy) + fabs(vz);
      
      vel[ 3*i+0 ] = vx;
      vel[ 3*i+1 ] = vy;
      vel[ 3*i+2 ] = vz;
#endif // BENCH
      
      lbmUpdateFlowField (i, density, vx, vy, vz, f_neq);
    }
  
  for (unit_level = 1; unit_level >= 0; unit_level--)
    {
      for (m = 0; m < net->neigh_procs; m++)
	{
	  neigh_proc_p = &net->neigh_proc[ m ];
	  
	  if ((unit_level == 1 && net->machine_id[ neigh_proc_p->id ] == net->machine_id[ net->id ]) ||
	      (unit_level == 0 && net->machine_id[ neigh_proc_p->id ] != net->machine_id[ net->id ]))
	    {
	      continue;
	    }
#ifndef NOMPI
	  net->err = MPI_Wait (&net->req[ 0 ][ m ], net->status);
	  net->err = MPI_Wait (&net->req[ 0 ][ net->neigh_procs + m ], net->status);
#endif
	  for (n = 0; n < neigh_proc_p->fs; n++)
	    {
	      f_new[ neigh_proc_p->f_recv_iv[n] ] = neigh_proc_p->f_to_recv[ n ];
	    }
	}
    }
  
  for (i = net->my_inner_sites; i < net->my_inner_sites + net->my_inter_sites; i++)
    {
      for (l = 0; l < 15; l++)
	{
#ifndef BENCH
	  if (f_old[ i*15+l ] < 0.) is_unstable = 1;
#endif // BENCH
	  
	  f_new[ f_id[i*15+l] ] = f_old[ i*15+l ];
	}
    }
  f_old_p = f_old;
  f_old = f_new;
  f_new = f_old_p;
  
#ifndef BENCH
  if (check_conv)
    {
      if (net->procs > 1)
	{
	  local_data[ 0 ] = 0.1 + (double)is_unstable;
	  local_data[ 1 ] = sum1;
	  local_data[ 2 ] = sum2;
#ifndef NOMPI
	  net->err = MPI_Allreduce (local_data, global_data, 3,
				    MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
	  global_data[ 0 ] = local_data[ 0 ];
	  global_data[ 1 ] = local_data[ 1 ];
	  global_data[ 2 ] = local_data[ 2 ];
#endif
	  is_unstable = (global_data[ 0 ] >= 1.);
	  sum1 = global_data[ 1 ];
	  sum2 = global_data[ 2 ];
	}
      
      if (sum1 <= sum2 * lbm->tolerance && sum2 > lbm->tolerance)
	{
	  *is_converged = 1;
	}
      lbm->conv_error = sum1 / sum2;
    }
  if (write_checkpoint)
    {
      lbmWriteConfig (!is_unstable, lbm->checkpoint_file_name, lbm, net);
    }
#endif // BENCH
  
  if (is_unstable)
    {
      return UNSTABLE;
    }
  else
    {
      return STABLE;
    }
}
*/

#ifndef TD
int lbmCycle (int write_checkpoint, int check_conv, int *is_converged, LBM *lbm, Net *net)
#else
int lbmCycle (int cycle_id, int time_step, int check_conv, int *is_converged, LBM *lbm, Net *net)
#endif
{
  // the entire simulation time step takes place through this function
  
  double f_neq[15];
  double omega;
  double density;
#ifndef TD
  double vx, vy, vz;
#else
  double vx[2], vy[2], vz[2];
#endif
  double *f_old_p;
  
#ifndef BENCH
  double sum1, sum2;
  double local_data[6];
  double global_data[6];
#endif // BENCH
  
#ifdef TD
  int is_current_cycle;
#endif
  int collision_type;
  int collision_count;
  int unit_level;
  int i, m, n;
  
  NeighProc *neigh_proc_p;
  
  
  for (m = 0; m < net->neigh_procs; m++)
    {
      neigh_proc_p = &net->neigh_proc[ m ];
      
#ifndef NOMPI
#ifndef TD
      n = neigh_proc_p->fs;
#else
      n = neigh_proc_p->fs * 2;
#endif
      net->err = MPI_Irecv (&neigh_proc_p->f_to_recv[ 0 ], n, MPI_DOUBLE,
			    neigh_proc_p->id, 10, MPI_COMM_WORLD,
			    &net->req[ 0 ][ net->neigh_procs + m ]);
#endif
    }
  
  *is_converged = 0;
  
#ifdef TD
  // If cycle_id is zero, only do the work for the current cycle.
  // Othewise, do the work for both cycles.
  if (cycle_id == 0)
    {
      is_current_cycle = 1;
      lbm->conv_error = 1.;
    }
  else
    {
      if (time_step == 0)
	{
	  lbm->conv_error = 0.;
	}
      is_current_cycle = 0;
    }
  if (time_step == 0)
    {
      lbm_density_min = +1.e+30;
      lbm_density_max = +1.e-30;
      
      lbm_velocity_min = +1.e+30;
      lbm_velocity_max = +1.e-30;
      
      lbm_stress_min = +1.e+30;
      lbm_stress_max = +1.e-30;
    }
#else
#ifndef BENCH
  lbm_density_min = +1.e+30;
  lbm_density_max = +1.e-30;
  
  lbm_velocity_min = +1.e+30;
  lbm_velocity_max = +1.e-30;
  
  lbm_stress_min = +1.e+30;
  lbm_stress_max = +1.e-30;
#endif
#endif
  
#ifndef BENCH
  is_unstable = 0;
  
  sum1 = 0.;
  sum2 = 0.;
#endif // BENCH
  
  omega = lbm->omega;
  
  collision_type = 0;
  collision_count = 0;
  
  for (i = net->my_inner_sites; i < net->my_sites; i++)
    {
#ifdef TD
      vx[0] = vy[0] = vz[0] = 1.e+30;
      
      for (is_current = is_current_cycle; is_current < 2; is_current++)
	{
	  (*lbmInterCollision[ collision_type ]) (omega, i, density,
						  &vx[is_current], &vy[is_current], &vz[is_current], f_neq);
	}
#ifndef BENCH
      sum1 += sqrt((vx[1] - vx[0]) * (vx[1] - vx[0]) +
		   (vy[1] - vy[0]) * (vy[1] - vy[0]) +
		   (vz[1] - vz[0]) * (vz[1] - vz[0]));
      sum2 += sqrt(vx[1] * vx[1] + vy[1] * vy[1] + vz[1] * vz[1]);
#endif // BENCH
#else // TD
      (*lbmInterCollision[ collision_type ]) (omega, i, &density, &vx, &vy, &vz, f_neq);
      
#ifndef BENCH
      sum1 += fabs(vel[3*i+0]-vx) + fabs(vel[3*i+1]-vy) + fabs(vel[3*i+2]-vz);
      sum2 += fabs(vx) + fabs(vy) + fabs(vz);
      
      vel[ 3*i+0 ] = vx;
      vel[ 3*i+1 ] = vy;
      vel[ 3*i+2 ] = vz;
#endif // BENCH
#endif // TD
      
#ifndef BENCH
#ifndef TD
      lbmUpdateFlowField (i, density, vx, vy, vz, f_neq);
#else
      lbmUpdateFlowField (i, density, vx[1], vy[1], vz[1], f_neq);
#endif
#endif // BENCH
      if (++collision_count == net->my_inter_collisions[ collision_type ])
	{
	  collision_count = 0;
	  
	  while (net->my_inter_collisions[ ++collision_type ] == 0) {;}
	}
    }
  
  for (unit_level = 1; unit_level >= 0; unit_level--)
    {
      for (m = 0; m < net->neigh_procs; m++)
	{
	  neigh_proc_p = &net->neigh_proc[ m ];
	  
	  if ((unit_level == 1 && net->machine_id[ neigh_proc_p->id ] == net->machine_id[ net->id ]) ||
	      (unit_level == 0 && net->machine_id[ neigh_proc_p->id ] != net->machine_id[ net->id ]))
	    {
	      continue;
	    }
	  for (n = 0; n < neigh_proc_p->fs; n++)
	    {
#ifndef TD
	      neigh_proc_p->f_to_send[ n ] = f_old[ neigh_proc_p->f_send_id[n] ];
#else
	      neigh_proc_p->f_to_send[ n*2   ] = f_old[ neigh_proc_p->f_send_id[n]    ];
	      neigh_proc_p->f_to_send[ n*2+1 ] = f_old[ neigh_proc_p->f_send_id[n]+15 ];
#endif
	    }
#ifndef NOMPI
#ifndef TD
	  n = neigh_proc_p->fs;
#else
	  n = neigh_proc_p->fs * 2;
#endif
	  net->err = MPI_Isend (&neigh_proc_p->f_to_send[ 0 ], n, MPI_DOUBLE,
				neigh_proc_p->id, 10, MPI_COMM_WORLD,
				&net->req[ 0 ][ m ]);
#endif
	}
    }

#ifdef NOOPENMP
  collision_type = 0;
  collision_count = 0;
#else
  int chunk = (int)ceil((double)net->my_inner_sites / (double)threads);
  
#ifndef BENCH
#pragma omp parallel for \
  default(shared) \
  reduction(+: sum1, sum2) \
  private(i, density, vx, vy, vz, f_neq, collision_type, collision_count) \
  schedule(dynamic, chunk)
#else // BENCH
#pragma omp parallel for \
  default(shared) \
  private(i, density, vx, vy, vz, f_neq, collision_type, collision_count) \
  schedule(dynamic, chunk)
#endif // BENCH
#endif // NOOPENMP
  for (i = 0; i < net->my_inner_sites; i++)
    {
#ifndef NOOPENMP
      collision_count = 0;
      
      for (collision_type = 0; collision_type < COLLISION_TYPES; collision_type++)
	{
	  collision_count += net->my_inner_collisions[ collision_type ];
	  
	  if (collision_count > i)
	    {
	      collision_count -= net->my_inner_collisions[ collision_type ];
	      collision_count = i - collision_count;
	      break;
	    }
	}
#endif // NOOPENMP

#ifdef TD
      vx[0] = vy[0] = vz[0] = 1.e+30;
      
      for (is_current = is_current_cycle; is_current < 2; is_current++)
	{
	  (*lbmInnerCollision[ collision_type ]) (omega, i, density,
						  vx[is_current], vy[is_current], vz[is_current], f_neq);
	}
#ifndef BENCH
      sum1 += sqrt((vx[1] - vx[0]) * (vx[1] - vx[0]) +
		   (vy[1] - vy[0]) * (vy[1] - vy[0]) +
		   (vz[1] - vz[0]) * (vz[1] - vz[0]));
      sum2 += sqrt(vx[1] * vx[1] + vy[1] * vy[1] + vz[1] * vz[1]);
#endif // BENCH
#else // TD
      (*lbmInnerCollision[ collision_type ]) (omega, i, &density, &vx, &vy, &vz, f_neq);
      
#ifndef BENCH
      sum1 += fabs(vel[ 3*i+0 ]-vx) + fabs(vel[ 3*i+1 ]-vy) + fabs(vel[ 3*i+2 ]-vz);
      sum2 += fabs(vx) + fabs(vy) + fabs(vz);
      
      vel[ 3*i+0 ] = vx;
      vel[ 3*i+1 ] = vy;
      vel[ 3*i+2 ] = vz;
#endif // BENCH
#endif // TD
      
#ifndef BENCH
#ifndef TD
      lbmUpdateFlowField (i, density, vx, vy, vz, f_neq);
#else
      lbmUpdateFlowField (i, density, vx[1], vy[1], vz[1], f_neq);
#endif
#endif // BENCH
      
#ifdef NOOPENMP
      if (++collision_count == net->my_inner_collisions[ collision_type ])
	{
	  collision_count = 0;
	  while (net->my_inner_collisions[ ++collision_type ] == 0) {;}
	}
#endif
    }
  
#ifndef BENCH
  if (check_conv)
    {
      if (net->procs > 1)
	{
	  local_data[ 0 ] = (double)is_unstable;
	  local_data[ 1 ] = sum1;
	  local_data[ 2 ] = sum2;
#ifndef NOMPI
	  net->err = MPI_Allreduce (local_data, global_data, 3,
				    MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
	  global_data[ 0 ] = local_data[ 0 ];
	  global_data[ 1 ] = local_data[ 1 ];
	  global_data[ 2 ] = local_data[ 2 ];
#endif
	  is_unstable = (global_data[ 0 ] >= 1.);
	  sum1 = global_data[ 1 ];
	  sum2 = global_data[ 2 ];
	}
    }
#ifndef TD
#ifndef NOMPI
  local_data[ 0 ] = lbm_density_min;
  local_data[ 1 ] = 1. / lbm_density_max;
  local_data[ 2 ] = lbm_velocity_min;
  local_data[ 3 ] = 1. / lbm_velocity_max;
  local_data[ 4 ] = lbm_stress_min;
  local_data[ 5 ] = 1. / lbm_stress_max;
  
  MPI_Reduce (local_data, global_data, 6,
  	      MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  
  lbm_density_min  = global_data[ 0 ];
  lbm_density_max  = 1. / global_data[ 1 ];
  lbm_velocity_min = global_data[ 2 ];
  lbm_velocity_max = 1. / global_data[ 3 ];
  lbm_stress_min   = global_data[ 4 ];
  lbm_stress_max   = 1. / global_data[ 5 ];
#endif
  if (check_conv)
    {
      if (sum1 <= sum2 * lbm->tolerance && sum2 > lbm->tolerance)
	{
	  *is_converged = 1;
	}
      lbm->conv_error = sum1 / sum2;
    }
  if (write_checkpoint)
    {
      lbmWriteConfig (!is_unstable, lbm->checkpoint_file_name, lbm, net);
    }
#else // TD
  if (time_step == lbm->period - 1)
    {
#ifndef NOMPI
      local_data[ 0 ] = lbm_density_min;
      local_data[ 1 ] = 1. / lbm_density_max;
      local_data[ 2 ] = lbm_velocity_min;
      local_data[ 3 ] = 1. / lbm_velocity_max;
      local_data[ 4 ] = lbm_stress_min;
      local_data[ 5 ] = 1. / lbm_stress_max;
      
      MPI_Reduce (local_data, global_data, 6,
		  MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
      
      lbm_density_min  = global_data[ 0 ];
      lbm_density_max  = 1. / global_data[ 1 ];
      lbm_velocity_min = global_data[ 2 ];
      lbm_velocity_max = 1. / global_data[ 3 ];
      lbm_stress_min   = global_data[ 4 ];
      lbm_stress_max   = 1. / global_data[ 5 ];
#endif
    }
  if (check_conv)
    {
      lbm->conv_error += sum1 / sum2;
      
      if (time_step == lbm->period - 1)
	{
	  lbm->conv_error /= lbm->period;
	  
	  if (lbm->conv_error < lbm->tolerance)
	    {
	      *is_converged = 1;
	    }
	}
    }
#endif // TD
#endif // BENCH
  
  for (unit_level = 1; unit_level >= 0; unit_level--)
    {
      for (m = 0; m < net->neigh_procs; m++)
	{
	  neigh_proc_p = &net->neigh_proc[ m ];
	  
	  if ((unit_level == 1 && net->machine_id[ neigh_proc_p->id ] == net->machine_id[ net->id ]) ||
	      (unit_level == 0 && net->machine_id[ neigh_proc_p->id ] != net->machine_id[ net->id ]))
	    {
	      continue;
	    }
#ifndef NOMPI
	  net->err = MPI_Wait (&net->req[ 0 ][ m ], net->status);
	  net->err = MPI_Wait (&net->req[ 0 ][ net->neigh_procs + m ], net->status);
#endif
	  for (n = 0; n < neigh_proc_p->fs; n++)
	    {
#ifndef TD
	      f_new[ neigh_proc_p->f_recv_iv[n] ] = neigh_proc_p->f_to_recv[ n ];
#else
	      f_new[ neigh_proc_p->f_recv_iv[n]    ] = neigh_proc_p->f_to_recv[ n*2   ];
	      f_new[ neigh_proc_p->f_recv_iv[n]+15 ] = neigh_proc_p->f_to_recv[ n*2+1 ];
#endif
	    }
	}
    }
  
  f_old_p = f_old;
  f_old = f_new;
  f_new = f_old_p;
  
  
#ifndef BENCH
  
  if (is_unstable)
    {
      return UNSTABLE;
    }
  else
    {
      return STABLE;
    }
#else // BENCH
  
  return STABLE;
#endif
}


void lbmEnd (LBM *lbm)
{
  free(outlet_density);
  outlet_density = NULL;
  
  free(inlet_density);
  inlet_density = NULL;
  
  free(lbm->fluid_sites_per_block);
  lbm->fluid_sites_per_block = NULL;
}
