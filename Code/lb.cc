// In this file, the functions useful to calculate the equilibrium distribution
// function, momentums, the effective von Mises stress and the boundary conditions
// are reported

// TODO: We shouldn't really be including config.h in here, but we have to at the moment.
#include "config.h"
#include "lb.h"
#include "utilityFunctions.h"

void (*lbmInnerCollision[COLLISION_TYPES]) (double omega, int i, double *density, double *v_x, double *v_y, double *v_z, double f_neq[]);
void (*lbmInterCollision[COLLISION_TYPES]) (double omega, int i, double *density, double *v_x, double *v_y, double *v_z, double f_neq[]);

void (*lbmUpdateSiteData[2][2]) (double omega, int i, double *density, double *vx, double *vy, double *vz, double *velocity,
				 void lbmCollision (double omega, int i,
						    double *density, double *v_x, double *v_y, double *v_z,
						    double f_neq[]));

double lbmConvertPressureToLatticeUnits (double pressure, LBM *lbm)
{
  return Cs2 + (pressure - REFERENCE_PRESSURE) * mmHg_TO_PASCAL *
    (PULSATILE_PERIOD / (lbm->period * lbm->voxel_size)) *
    (PULSATILE_PERIOD / (lbm->period * lbm->voxel_size)) / BLOOD_DENSITY;
}


double lbmConvertPressureToPhysicalUnits (double pressure, LBM *lbm)
{
  return REFERENCE_PRESSURE + ((pressure / Cs2 - 1.0) * Cs2) * BLOOD_DENSITY *
    ((lbm->period * lbm->voxel_size) / PULSATILE_PERIOD) *
    ((lbm->period * lbm->voxel_size) / PULSATILE_PERIOD) / mmHg_TO_PASCAL;
}


double lbmConvertPressureGradToLatticeUnits (double pressure_grad, LBM *lbm)
{
  return pressure_grad * mmHg_TO_PASCAL *
    (PULSATILE_PERIOD / (lbm->period * lbm->voxel_size)) *
    (PULSATILE_PERIOD / (lbm->period * lbm->voxel_size)) / BLOOD_DENSITY;
}


double lbmConvertPressureGradToPhysicalUnits (double pressure_grad, LBM *lbm)
{
  return pressure_grad * BLOOD_DENSITY *
    ((lbm->period * lbm->voxel_size) / PULSATILE_PERIOD) *
    ((lbm->period * lbm->voxel_size) / PULSATILE_PERIOD) / mmHg_TO_PASCAL;
}


double lbmConvertVelocityToLatticeUnits (double velocity, LBM *lbm)
{
  return velocity * (((lbm->tau-0.5)/3.0) * lbm->voxel_size) / (BLOOD_VISCOSITY / BLOOD_DENSITY);
}


double lbmConvertVelocityToPhysicalUnits (double velocity, LBM *lbm)
{
  // convert velocity from lattice units to physical units (m/s)
  
  return velocity * (BLOOD_VISCOSITY / BLOOD_DENSITY) / (((lbm->tau-0.5)/3.0) * lbm->voxel_size);
}


double lbmConvertStressToLatticeUnits (double stress, LBM *lbm)
{
  return stress * (BLOOD_DENSITY / (BLOOD_VISCOSITY * BLOOD_VISCOSITY)) *
    (((lbm->tau-0.5)/3.0) * lbm->voxel_size) *
    (((lbm->tau-0.5)/3.0) * lbm->voxel_size);
}


double lbmConvertStressToPhysicalUnits (double stress, LBM *lbm)
{
  // convert stress from lattice units to physical units (Pa)
  
  return stress * BLOOD_VISCOSITY * BLOOD_VISCOSITY /
    (BLOOD_DENSITY * (((lbm->tau-0.5)/3.0) * lbm->voxel_size) * (((lbm->tau-0.5)/3.0) * lbm->voxel_size));
}


double lbmCalculateTau (LBM *lbm)
{
  return 0.5 + (PULSATILE_PERIOD * BLOOD_VISCOSITY / BLOOD_DENSITY) /
    (Cs2 * lbm->period * lbm->voxel_size * lbm->voxel_size);
}


// Calculate density, velocity and the equilibrium distribution
// functions according to the D3Q15 model.  The calculated v_x, v_y
// and v_z are actually density * velocity, because we are using the
// compressible model.
void lbmFeq (double f[], double *density, double *v_x, double *v_y, double *v_z, double f_eq[])
{
  double density_1;
  double v_xx, v_yy, v_zz;
  double temp1, temp2;
  
  
  *v_x = f[1] + (f[7] + f[9]) + (f[11] + f[13]);
  *v_y = f[3] + (f[12] + f[14]);
  *v_z = f[5] + f[10];
  
  *density = f[0] + (f[2] + f[4]) + (f[6] + f[8]) + *v_x + *v_y + *v_z;
  
  *v_x -= f[2] + (f[8] + f[10]) + (f[12] + f[14]);
  *v_y += (f[7] + f[9]) - (f[4] + (f[8] + f[10]) + (f[11] + f[13]));
  *v_z += f[7] + f[11] + f[14] - ((f[6] + f[8]) + f[9] + f[12] + f[13]);
  
  density_1 = 1. / *density;
  
  v_xx = *v_x * *v_x;
  v_yy = *v_y * *v_y;
  v_zz = *v_z * *v_z;
  
  f_eq[ 0 ] = (2.0/9.0) * *density - (1.0/3.0) * ((v_xx + v_yy + v_zz) * density_1);
  
  temp1 = (1.0/9.0) * *density - (1.0/6.0) * ((v_xx + v_yy + v_zz) * density_1);
  
  f_eq[1] = (temp1 + (0.5 * density_1) * v_xx) + ((1.0/3.0) * *v_x);   // (+1, 0, 0)
  f_eq[2] = (temp1 + (0.5 * density_1) * v_xx) - ((1.0/3.0) * *v_x);   // (+1, 0, 0)
  
  f_eq[3] = (temp1 + (0.5 * density_1) * v_yy) + ((1.0/3.0) * *v_y);   // (0, +1, 0)
  f_eq[4] = (temp1 + (0.5 * density_1) * v_yy) - ((1.0/3.0) * *v_y);   // (0, -1, 0)
  
  f_eq[5] = (temp1 + (0.5 * density_1) * v_zz) + ((1.0/3.0) * *v_z);   // (0, 0, +1)
  f_eq[6] = (temp1 + (0.5 * density_1) * v_zz) - ((1.0/3.0) * *v_z);   // (0, 0, -1)
  
  temp1 *= (1.0/8.0);
  
  temp2 = (*v_x + *v_y) + *v_z;
  
  f_eq[ 7] = (temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) + ((1.0/24.0) * temp2);   // (+1, +1, +1)
  f_eq[ 8] = (temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) - ((1.0/24.0) * temp2);   // (-1, -1, -1)
  
  temp2 = (*v_x + *v_y) - *v_z;
  
  f_eq[ 9] = (temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) + ((1.0/24.0) * temp2);   // (+1, +1, -1)
  f_eq[10] = (temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) - ((1.0/24.0) * temp2);   // (-1, -1, +1)
  
  temp2 = (*v_x - *v_y) + *v_z;
  
  f_eq[11] = (temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) + ((1.0/24.0) * temp2);   // (+1, -1, +1)
  f_eq[12] = (temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) - ((1.0/24.0) * temp2);   // (-1, +1, -1)
  
  temp2 = (*v_x - *v_y) - *v_z;
  
  f_eq[13] = (temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) + ((1.0/24.0) * temp2);   // (+1, -1, -1)
  f_eq[14] = (temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) - ((1.0/24.0) * temp2);   // (-1, +1, +1)
}


void lbmFeq (double density, double v_x, double v_y, double v_z, double f_eq[])
{
  double density_1;
  double v_xx, v_yy, v_zz;
  double temp1, temp2;
  
  
  density_1 = 1. / density;
  
  v_xx = v_x * v_x;
  v_yy = v_y * v_y;
  v_zz = v_z * v_z;
  
  f_eq[ 0 ] = (2.0/9.0) * density - (1.0/3.0) * ((v_xx + v_yy + v_zz) * density_1);
  
  temp1 = (1.0/9.0) * density - (1.0/6.0) * ((v_xx + v_yy + v_zz) * density_1);
  
  f_eq[1] = (temp1 + (0.5 * density_1) * v_xx) + ((1.0/3.0) * v_x);   // (+1, 0, 0)
  f_eq[2] = (temp1 + (0.5 * density_1) * v_xx) - ((1.0/3.0) * v_x);   // (+1, 0, 0)
  
  f_eq[3] = (temp1 + (0.5 * density_1) * v_yy) + ((1.0/3.0) * v_y);   // (0, +1, 0)
  f_eq[4] = (temp1 + (0.5 * density_1) * v_yy) - ((1.0/3.0) * v_y);   // (0, -1, 0)
  
  f_eq[5] = (temp1 + (0.5 * density_1) * v_zz) + ((1.0/3.0) * v_z);   // (0, 0, +1)
  f_eq[6] = (temp1 + (0.5 * density_1) * v_zz) - ((1.0/3.0) * v_z);   // (0, 0, -1)
  
  temp1 *= (1.0/8.0);
  
  temp2 = (v_x + v_y) + v_z;
  
  f_eq[ 7] = (temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) + ((1.0/24.0) * temp2);   // (+1, +1, +1)
  f_eq[ 8] = (temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) - ((1.0/24.0) * temp2);   // (-1, -1, -1)
  
  temp2 = (v_x + v_y) - v_z;
  
  f_eq[ 9] = (temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) + ((1.0/24.0) * temp2);   // (+1, +1, -1)
  f_eq[10] = (temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) - ((1.0/24.0) * temp2);   // (-1, -1, +1)
  
  temp2 = (v_x - v_y) + v_z;
  
  f_eq[11] = (temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) + ((1.0/24.0) * temp2);   // (+1, -1, +1)
  f_eq[12] = (temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) - ((1.0/24.0) * temp2);   // (-1, +1, -1)
  
  temp2 = (v_x - v_y) - v_z;
  
  f_eq[13] = (temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) + ((1.0/24.0) * temp2);   // (+1, -1, -1)
  f_eq[14] = (temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) - ((1.0/24.0) * temp2);   // (-1, +1, +1)
}


// Collision + streaming for non-boundary fluid lattice sites non-adjacent
// to neigbhouring subdomains.
void lbmInnerCollision0 (double omega, int i,
			 double *density, double *v_x, double *v_y, double *v_z,
			 double f_neq[])
{
  double density_1;
  double *f;
  double v_xx, v_yy, v_zz;
  double temp1, temp2;
  
  
  f = &f_old[ i*15 ];
  
  *v_x = f[1] + (f[7] + f[9]) + (f[11] + f[13]);
  *v_y = f[3] + (f[12] + f[14]);
  *v_z = f[5] + f[10];
  
  *density = f[0] + (f[2] + f[4]) + (f[6] + f[ 8 ]) + *v_x + *v_y + *v_z;
  
  *v_x -= f[2] + (f[8] + f[10]) + (f[12] + f[14]);
  *v_y += (f[7] + f[9]) - (f[4] + (f[8] + f[10]) + (f[11] + f[13]));
  *v_z += f[7] + f[11] + f[14] - ((f[6] + f[8]) + f[9] + f[12] + f[13]);
  
  density_1 = 1. / *density;
  
  v_xx = *v_x * *v_x;
  v_yy = *v_y * *v_y;
  v_zz = *v_z * *v_z;
  
  f_new[ f_id[i*15+0] ] = f[0] + omega * (f_neq[0] = f[0] - ((2.0/9.0) * *density - (1.0/3.0) * ((v_xx + v_yy + v_zz) * density_1)));
  
  temp1 = (1.0/9.0) * *density - (1.0/6.0) * ((v_xx + v_yy + v_zz) * density_1);
  
  f_new[ f_id[i*15+1] ] = f[1] + omega * (f_neq[1] = f[1] - ((temp1 + (0.5 * density_1) * v_xx) + ((1.0/3.0) * *v_x)));   // (+1, 0, 0)
  f_new[ f_id[i*15+2] ] = f[2] + omega * (f_neq[2] = f[2] - ((temp1 + (0.5 * density_1) * v_xx) - ((1.0/3.0) * *v_x)));   // (+1, 0, 0)
  
  f_new[ f_id[i*15+3] ] = f[3] + omega * (f_neq[3] = f[3] - ((temp1 + (0.5 * density_1) * v_yy) + ((1.0/3.0) * *v_y)));   // (0, +1, 0)
  f_new[ f_id[i*15+4] ] = f[4] + omega * (f_neq[4] = f[4] - ((temp1 + (0.5 * density_1) * v_yy) - ((1.0/3.0) * *v_y)));   // (0, -1, 0)
  
  f_new[ f_id[i*15+5] ] = f[5] + omega * (f_neq[5] = f[5] - ((temp1 + (0.5 * density_1) * v_zz) + ((1.0/3.0) * *v_z)));   // (0, 0, +1)
  f_new[ f_id[i*15+6] ] = f[6] + omega * (f_neq[6] = f[6] - ((temp1 + (0.5 * density_1) * v_zz) - ((1.0/3.0) * *v_z)));   // (0, 0, -1)
  
  temp1 *= (1.0/8.0);
  
  temp2 = (*v_x + *v_y) + *v_z;
  
  f_new[ f_id[i*15+7] ] = f[7] + omega * (f_neq[7] = f[7] - ((temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) + ((1.0/24.0) * temp2)));   // (+1, +1, +1)
  f_new[ f_id[i*15+8] ] = f[8] + omega * (f_neq[8] = f[8] - ((temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) - ((1.0/24.0) * temp2)));   // (-1, -1, -1)
  
  temp2 = (*v_x + *v_y) - *v_z;
  
  f_new[ f_id[i*15+ 9] ] = f[ 9] + omega * (f_neq[ 9] = f[ 9] - ((temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) + ((1.0/24.0) * temp2)));   // (+1, +1, -1)
  f_new[ f_id[i*15+10] ] = f[10] + omega * (f_neq[10] = f[10] - ((temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) - ((1.0/24.0) * temp2)));   // (-1, -1, +1)
  
  temp2 = (*v_x - *v_y) + *v_z;
  
  f_new[ f_id[i*15+11] ] = f[11] + omega * (f_neq[11] = f[11] - ((temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) + ((1.0/24.0) * temp2)));   // (+1, -1, +1)
  f_new[ f_id[i*15+12] ] = f[12] + omega * (f_neq[12] = f[12] - ((temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) - ((1.0/24.0) * temp2)));   // (-1, +1, -1)
  
  temp2 = (*v_x - *v_y) - *v_z;	 
  
  f_new[ f_id[i*15+13] ] = f[13] + omega * (f_neq[13] = f[13] - ((temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) + ((1.0/24.0) * temp2)));   // (+1, -1, -1)
  f_new[ f_id[i*15+14] ] = f[14] + omega * (f_neq[14] = f[14] - ((temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) - ((1.0/24.0) * temp2)));   // (-1, +1, +1)
}


// Collision + streaming for non-boundary fluid lattice sites adjacent
// to neigbhouring subdomains.
void lbmInterCollision0 (double omega, int i,
			 double *density, double *v_x, double *v_y, double *v_z,
			 double f_neq[])
{
  double density_1;
  double *f;
  double v_xx, v_yy, v_zz;
  double temp1, temp2;
  
  
  f = &f_old[ i*15 ];
  
  *v_x = f[1] + (f[7] + f[9]) + (f[11] + f[13]);
  *v_y = f[3] + (f[12] + f[14]);
  *v_z = f[5] + f[10];
  
  *density = f[0] + (f[2] + f[4]) + (f[6] + f[ 8 ]) + *v_x + *v_y + *v_z;
  
  *v_x -= f[2] + (f[8] + f[10]) + (f[12] + f[14]);
  *v_y += (f[7] + f[9]) - (f[4] + (f[8] + f[10]) + (f[11] + f[13]));
  *v_z += f[7] + f[11] + f[14] - ((f[6] + f[8]) + f[9] + f[12] + f[13]);
  
  density_1 = 1. / *density;
  
  v_xx = *v_x * *v_x;
  v_yy = *v_y * *v_y;
  v_zz = *v_z * *v_z;
  
  f_new[ f_id[i*15+0] ] = f[0] + omega * (f_neq[0] = f[0] - ((2.0/9.0) * *density - (1.0/3.0) * ((v_xx + v_yy + v_zz) * density_1)));
  
  temp1 = (1.0/9.0) * *density - (1.0/6.0) * ((v_xx + v_yy + v_zz) * density_1);
  
  f_new[ f_id[i*15+1] ] = f[1] += omega * (f_neq[1] = f[1] - ((temp1 + (0.5 * density_1) * v_xx) + ((1.0/3.0) * *v_x)));   // (+1, 0, 0)
  f_new[ f_id[i*15+2] ] = f[2] += omega * (f_neq[2] = f[2] - ((temp1 + (0.5 * density_1) * v_xx) - ((1.0/3.0) * *v_x)));   // (+1, 0, 0)
  
  f_new[ f_id[i*15+3] ] = f[3] += omega * (f_neq[3] = f[3] - ((temp1 + (0.5 * density_1) * v_yy) + ((1.0/3.0) * *v_y)));   // (0, +1, 0)
  f_new[ f_id[i*15+4] ] = f[4] += omega * (f_neq[4] = f[4] - ((temp1 + (0.5 * density_1) * v_yy) - ((1.0/3.0) * *v_y)));   // (0, -1, 0)
  
  f_new[ f_id[i*15+5] ] = f[5] += omega * (f_neq[5] = f[5] - ((temp1 + (0.5 * density_1) * v_zz) + ((1.0/3.0) * *v_z)));   // (0, 0, +1)
  f_new[ f_id[i*15+6] ] = f[6] += omega * (f_neq[6] = f[6] - ((temp1 + (0.5 * density_1) * v_zz) - ((1.0/3.0) * *v_z)));   // (0, 0, -1)
  
  temp1 *= (1.0/8.0);
  
  temp2 = (*v_x + *v_y) + *v_z;
  
  f_new[ f_id[i*15+7] ] = f[7] += omega * (f_neq[7] = f[7] - ((temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) + ((1.0/24.0) * temp2)));   // (+1, +1, +1)
  f_new[ f_id[i*15+8] ] = f[8] += omega * (f_neq[8] = f[8] - ((temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) - ((1.0/24.0) * temp2)));   // (-1, -1, -1)
  
  temp2 = (*v_x + *v_y) - *v_z;
  
  f_new[ f_id[i*15+ 9] ] = f[ 9] += omega * (f_neq[ 9] = f[ 9] - ((temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) + ((1.0/24.0) * temp2)));   // (+1, +1, -1)
  f_new[ f_id[i*15+10] ] = f[10] += omega * (f_neq[10] = f[10] - ((temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) - ((1.0/24.0) * temp2)));   // (-1, -1, +1)
  
  temp2 = (*v_x - *v_y) + *v_z;
  
  f_new[ f_id[i*15+11] ] = f[11] += omega * (f_neq[11] = f[11] - ((temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) + ((1.0/24.0) * temp2)));   // (+1, -1, +1)
  f_new[ f_id[i*15+12] ] = f[12] += omega * (f_neq[12] = f[12] - ((temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) - ((1.0/24.0) * temp2)));   // (-1, +1, -1)
  
  temp2 = (*v_x - *v_y) - *v_z;	 
  
  f_new[ f_id[i*15+13] ] = f[13] += omega * (f_neq[13] = f[13] - ((temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) + ((1.0/24.0) * temp2)));   // (+1, -1, -1)
  f_new[ f_id[i*15+14] ] = f[14] += omega * (f_neq[14] = f[14] - ((temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) - ((1.0/24.0) * temp2)));   // (-1, +1, +1)
}


// Collision + streaming for fluid lattice sites adjacent to the wall.
void lbmCollision1 (double omega, int i,
		    double *density, double *v_x, double *v_y, double *v_z,
		    double f_neq[])
{
  double *f;
  double temp;
  
  int l;
  
  
  f = &f_old[ i*15 ];
  
  for (l = 0; l < 15; l++)
    {
      f_neq[l] = f[l];
    }
  *v_x = *v_y = *v_z = 0.F;
  
  *density = 0.;
  
  for (l = 0; l < 15; l++) *density += f[l];
  
  f_neq[0] -= (f_new[ f_id[i*15] ] = f[0] = (2.0/9.0) * *density);
  
  temp = (1.0/9.0) * *density;
  
  for (l = 1; l < 7; l++)
    f_neq[l] -= (f_new[ f_id[i*15+l] ] = f[l] = temp);
  
  temp *= (1.0/8.0);
  
  for (l = 7; l < 15; l++)
    f_neq[l] -= (f_new[ f_id[i*15+l] ] = f[l] = temp);
}


// Collision + streaming for inlet fluid lattice sites.
void lbmCollision2 (double omega, int i,
		    double *density, double *v_x, double *v_y, double *v_z,
		    double f_neq[])
{
  double *f;
  double dummy_density;
  
  unsigned int boundary_id, l;
  
  
  f = &f_old[ i*15 ];
  
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
      f_neq[l] -= (f_new[ f_id[i*15+l] ] = f[l]);
    }
}


// Collision + streaming for outlet fluid lattice sites.
void lbmCollision3 (double omega, int i,
		    double *density, double *v_x, double *v_y, double *v_z,
		    double f_neq[])
{
  double *f;
  double dummy_density;
  
  unsigned int boundary_id, l;
  
  
  f = &f_old[ i*15 ];
  
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
      f_neq[l] -= (f_new[ f_id[i*15+l] ] = f[l]);
    }
}


// Collision + streaming for fluid lattice sites and adjacent to the inlet and the wall.
void lbmCollision4 (double omega, int i,
		    double *density, double *v_x, double *v_y, double *v_z,
		    double f_neq[])
{
  double *f;
  double temp;
  
  int l;
  
  unsigned int boundary_id;
  
  
  f = &f_old[ i*15 ];
  
  for (l = 0; l < 15; l++)
    {
      f_neq[l] = f[l];
    }
  boundary_id = (net_site_data[ i ] & BOUNDARY_ID_MASK) >> BOUNDARY_ID_SHIFT;
  
  *density = inlet_density[ boundary_id ];
  
  *v_x = *v_y = *v_z = 0.F;
  
  f_neq[0] -= (f_new[ f_id[i*15] ] = f[0] = (2.0/9.0) * *density);
  
  temp = (1.0/9.0) * *density;
  
  for (l = 1; l < 7; l++)
    f_neq[l] -= (f_new[ f_id[i*15+l] ] = f[l] = temp);
  
  temp *= (1.0/8.0);
  
  for (l = 7; l < 15; l++)
    f_neq[l] -= (f_new[ f_id[i*15+l] ] = f[l] = temp);
}


// Collision + streaming for fluid lattice sites and adjacent to the outlet and the wall.
void lbmCollision5 (double omega, int i,
		    double *density, double *v_x, double *v_y, double *v_z,
		    double f_neq[])
{
  double *f;
  double temp;
  
  int l;
  
  unsigned int boundary_id;
  
  
  f = &f_old[ i*15 ];
  
  for (l = 0; l < 15; l++)
    {
      f_neq[l] = f[l];
    }
  boundary_id = (net_site_data[ i ] & BOUNDARY_ID_MASK) >> BOUNDARY_ID_SHIFT;
  
  *density = outlet_density[ boundary_id ];
  
  *v_x = *v_y = *v_z = 0.F;
  
  f_neq[0] -= (f_new[ f_id[i*15] ] = f[0] = (2.0/9.0) * *density);
  
  temp = (1.0/9.0) * *density;
  
  for (l = 1; l < 7; l++)
    f_neq[l] -= (f_new[ f_id[i*15+l] ] = f[l] = temp);
  
  temp *= (1.0/8.0);
  
  for (l = 7; l < 15; l++)
    f_neq[l] -= (f_new[ f_id[i*15+l] ] = f[l] = temp);
}


// The same as lbmInnerCollision0 but useful for convergence purposes.
void lbmInnerCollisionConv0 (double omega, int i,
			     double *density, double *v_x, double *v_y, double *v_z,
			     double f_neq[])
{
  double density_1;
  double *f;
  double v_xx, v_yy, v_zz;
  double temp1, temp2;
  
  
  f = &f_old[ i*30+cycle_tag*15 ];
  
  *v_x = f[1] + (f[7] + f[9]) + (f[11] + f[13]);
  *v_y = f[3] + (f[12] + f[14]);
  *v_z = f[5] + f[10];
  
  *density = f[0] + (f[2] + f[4]) + (f[6] + f[ 8 ]) + *v_x + *v_y + *v_z;
  
  *v_x -= f[2] + (f[8] + f[10]) + (f[12] + f[14]);
  *v_y += (f[7] + f[9]) - (f[4] + (f[8] + f[10]) + (f[11] + f[13]));
  *v_z += f[7] + f[11] + f[14] - ((f[6] + f[8]) + f[9] + f[12] + f[13]);
  
  density_1 = 1. / *density;
  
  v_xx = *v_x * *v_x;
  v_yy = *v_y * *v_y;
  v_zz = *v_z * *v_z;
  
  f_new[ f_id[i*15+0]+cycle_tag*15 ] = f[0] + omega * (f_neq[0] = f[0] - ((2.0/9.0) * *density - (1.0/3.0) * ((v_xx + v_yy + v_zz) * density_1)));
  
  temp1 = (1.0/9.0) * *density - (1.0/6.0) * ((v_xx + v_yy + v_zz) * density_1);
  
  f_new[ f_id[i*15+1]+cycle_tag*15 ] = f[1] + omega * (f_neq[1] = f[1] - ((temp1 + (0.5 * density_1) * v_xx) + ((1.0/3.0) * *v_x)));   // (+1, 0, 0)
  f_new[ f_id[i*15+2]+cycle_tag*15 ] = f[2] + omega * (f_neq[2] = f[2] - ((temp1 + (0.5 * density_1) * v_xx) - ((1.0/3.0) * *v_x)));   // (+1, 0, 0)
  
  f_new[ f_id[i*15+3]+cycle_tag*15 ] = f[3] + omega * (f_neq[3] = f[3] - ((temp1 + (0.5 * density_1) * v_yy) + ((1.0/3.0) * *v_y)));   // (0, +1, 0)
  f_new[ f_id[i*15+4]+cycle_tag*15 ] = f[4] + omega * (f_neq[4] = f[4] - ((temp1 + (0.5 * density_1) * v_yy) - ((1.0/3.0) * *v_y)));   // (0, -1, 0)
  
  f_new[ f_id[i*15+5]+cycle_tag*15 ] = f[5] + omega * (f_neq[5] = f[5] - ((temp1 + (0.5 * density_1) * v_zz) + ((1.0/3.0) * *v_z)));   // (0, 0, +1)
  f_new[ f_id[i*15+6]+cycle_tag*15 ] = f[6] + omega * (f_neq[6] = f[6] - ((temp1 + (0.5 * density_1) * v_zz) - ((1.0/3.0) * *v_z)));   // (0, 0, -1)
  
  temp1 *= (1.0/8.0);
  
  temp2 = (*v_x + *v_y) + *v_z;
  
  f_new[ f_id[i*15+7]+cycle_tag*15 ] = f[7] + omega * (f_neq[7] = f[7] - ((temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) + ((1.0/24.0) * temp2)));   // (+1, +1, +1)
  f_new[ f_id[i*15+8]+cycle_tag*15 ] = f[8] + omega * (f_neq[8] = f[8] - ((temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) - ((1.0/24.0) * temp2)));   // (-1, -1, -1)
  
  temp2 = (*v_x + *v_y) - *v_z;
  
  f_new[ f_id[i*15+ 9]+cycle_tag*15 ] = f[ 9] + omega * (f_neq[ 9] = f[ 9] - ((temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) + ((1.0/24.0) * temp2)));   // (+1, +1, -1)
  f_new[ f_id[i*15+10]+cycle_tag*15 ] = f[10] + omega * (f_neq[10] = f[10] - ((temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) - ((1.0/24.0) * temp2)));   // (-1, -1, +1)
  
  temp2 = (*v_x - *v_y) + *v_z;
  
  f_new[ f_id[i*15+11]+cycle_tag*15 ] = f[11] + omega * (f_neq[11] = f[11] - ((temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) + ((1.0/24.0) * temp2)));   // (+1, -1, +1)
  f_new[ f_id[i*15+12]+cycle_tag*15 ] = f[12] + omega * (f_neq[12] = f[12] - ((temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) - ((1.0/24.0) * temp2)));   // (-1, +1, -1)
  
  temp2 = (*v_x - *v_y) - *v_z;	 
  
  f_new[ f_id[i*15+13]+cycle_tag*15 ] = f[13] + omega * (f_neq[13] = f[13] - ((temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) + ((1.0/24.0) * temp2)));   // (+1, -1, -1)
  f_new[ f_id[i*15+14]+cycle_tag*15 ] = f[14] + omega * (f_neq[14] = f[14] - ((temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) - ((1.0/24.0) * temp2)));   // (-1, +1, +1)
}


// The same as lbmInterCollision0 but useful for convergence purposes.
void lbmInterCollisionConv0 (double omega, int i,
			     double *density, double *v_x, double *v_y, double *v_z,
			     double f_neq[])
{
  double density_1;
  double *f;
  double v_xx, v_yy, v_zz;
  double temp1, temp2;
  
  
  f = &f_old[ i*30+cycle_tag*15 ];
  
  *v_x = f[1] + (f[7] + f[9]) + (f[11] + f[13]);
  *v_y = f[3] + (f[12] + f[14]);
  *v_z = f[5] + f[10];
  
  *density = f[0] + (f[2] + f[4]) + (f[6] + f[ 8 ]) + *v_x + *v_y + *v_z;
  
  *v_x -= f[2] + (f[8] + f[10]) + (f[12] + f[14]);
  *v_y += (f[7] + f[9]) - (f[4] + (f[8] + f[10]) + (f[11] + f[13]));
  *v_z += f[7] + f[11] + f[14] - ((f[6] + f[8]) + f[9] + f[12] + f[13]);
  
  density_1 = 1. / *density;
  
  v_xx = *v_x * *v_x;
  v_yy = *v_y * *v_y;
  v_zz = *v_z * *v_z;
  
  f_new[ f_id[i*15+0]+cycle_tag*15 ] = f[0] + omega * (f_neq[0] = f[0] - ((2.0/9.0) * *density - (1.0/3.0) * ((v_xx + v_yy + v_zz) * density_1)));
  
  temp1 = (1.0/9.0) * *density - (1.0/6.0) * ((v_xx + v_yy + v_zz) * density_1);
  
  f_new[ f_id[i*15+1]+cycle_tag*15 ] = f[1] += omega * (f_neq[1] = f[1] - ((temp1 + (0.5 * density_1) * v_xx) + ((1.0/3.0) * *v_x)));   // (+1, 0, 0)
  f_new[ f_id[i*15+2]+cycle_tag*15 ] = f[2] += omega * (f_neq[2] = f[2] - ((temp1 + (0.5 * density_1) * v_xx) - ((1.0/3.0) * *v_x)));   // (+1, 0, 0)
  
  f_new[ f_id[i*15+3]+cycle_tag*15 ] = f[3] += omega * (f_neq[3] = f[3] - ((temp1 + (0.5 * density_1) * v_yy) + ((1.0/3.0) * *v_y)));   // (0, +1, 0)
  f_new[ f_id[i*15+4]+cycle_tag*15 ] = f[4] += omega * (f_neq[4] = f[4] - ((temp1 + (0.5 * density_1) * v_yy) - ((1.0/3.0) * *v_y)));   // (0, -1, 0)
  
  f_new[ f_id[i*15+5]+cycle_tag*15 ] = f[5] += omega * (f_neq[5] = f[5] - ((temp1 + (0.5 * density_1) * v_zz) + ((1.0/3.0) * *v_z)));   // (0, 0, +1)
  f_new[ f_id[i*15+6]+cycle_tag*15 ] = f[6] += omega * (f_neq[6] = f[6] - ((temp1 + (0.5 * density_1) * v_zz) - ((1.0/3.0) * *v_z)));   // (0, 0, -1)
  
  temp1 *= (1.0/8.0);
  
  temp2 = (*v_x + *v_y) + *v_z;
  
  f_new[ f_id[i*15+7]+cycle_tag*15 ] = f[7] += omega * (f_neq[7] = f[7] - ((temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) + ((1.0/24.0) * temp2)));   // (+1, +1, +1)
  f_new[ f_id[i*15+8]+cycle_tag*15 ] = f[8] += omega * (f_neq[8] = f[8] - ((temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) - ((1.0/24.0) * temp2)));   // (-1, -1, -1)
  
  temp2 = (*v_x + *v_y) - *v_z;
  
  f_new[ f_id[i*15+ 9]+cycle_tag*15 ] = f[ 9] += omega * (f_neq[ 9] = f[ 9] - ((temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) + ((1.0/24.0) * temp2)));   // (+1, +1, -1)
  f_new[ f_id[i*15+10]+cycle_tag*15 ] = f[10] += omega * (f_neq[10] = f[10] - ((temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) - ((1.0/24.0) * temp2)));   // (-1, -1, +1)
  
  temp2 = (*v_x - *v_y) + *v_z;
  
  f_new[ f_id[i*15+11]+cycle_tag*15 ] = f[11] += omega * (f_neq[11] = f[11] - ((temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) + ((1.0/24.0) * temp2)));   // (+1, -1, +1)
  f_new[ f_id[i*15+12]+cycle_tag*15 ] = f[12] += omega * (f_neq[12] = f[12] - ((temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) - ((1.0/24.0) * temp2)));   // (-1, +1, -1)
  
  temp2 = (*v_x - *v_y) - *v_z;	 
  
  f_new[ f_id[i*15+13]+cycle_tag*15 ] = f[13] += omega * (f_neq[13] = f[13] - ((temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) + ((1.0/24.0) * temp2)));   // (+1, -1, -1)
  f_new[ f_id[i*15+14]+cycle_tag*15 ] = f[14] += omega * (f_neq[14] = f[14] - ((temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) - ((1.0/24.0) * temp2)));   // (-1, +1, +1)
}


// The same as lbmCollision1 but useful for convergence purposes.
void lbmCollisionConv1 (double omega, int i,
			double *density, double *v_x, double *v_y, double *v_z,
			double f_neq[])
{
  double *f;
  double temp;
  
  int l;
  
  
  f = &f_old[ i*30+cycle_tag*15 ];
  
  for (l = 0; l < 15; l++)
    {
      f_neq[l] = f[l];
    }
  *v_x = *v_y = *v_z = 0.F;
  
  *density = 0.;
  
  for (l = 0; l < 15; l++) *density += f[l];
  
  f_neq[0] -= (f_new[ f_id[i*15]+cycle_tag*15 ] = f[0] = (2.0/9.0) * *density);
  
  temp = (1.0/9.0) * *density;
  
  for (l = 1; l < 7; l++)
    f_neq[l] -= (f_new[ f_id[i*15+l]+cycle_tag*15 ] = f[l] = temp);
  
  temp *= (1.0/8.0);
  
  for (l = 7; l < 15; l++)
    f_neq[l] -= (f_new[ f_id[i*15+l]+cycle_tag*15 ] = f[l] = temp);
}


// The same as lbmCollision2 but useful for convergence purposes.
void lbmCollisionConv2 (double omega, int i,
			double *density, double *v_x, double *v_y, double *v_z,
			double f_neq[])
{
  double *f;
  double dummy_density;
  
  unsigned int boundary_id, l;
  
  
  f = &f_old[ i*30+cycle_tag*15 ];
  
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
      f_neq[l] -= (f_new[ f_id[i*15+l]+cycle_tag*15 ] = f[l]);
    }
}


// The same as lbmCollision3 but useful for convergence purposes.
void lbmCollisionConv3 (double omega, int i,
			double *density, double *v_x, double *v_y, double *v_z,
			double f_neq[])
{
  double *f;
  double dummy_density;
  
  unsigned int boundary_id, l;
  
  
  f = &f_old[ i*30+cycle_tag*15 ];
  
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
      f_neq[l] -= (f_new[ f_id[i*15+l]+cycle_tag*15 ] = f[l]);
    }
}


// The same as lbmCollision4 but useful for convergence purposes.
void lbmCollisionConv4 (double omega, int i,
			double *density, double *v_x, double *v_y, double *v_z,
			double f_neq[])
{
  double *f;
  double temp;
  
  int l;
  
  unsigned int boundary_id;
  
  
  f = &f_old[ i*30+cycle_tag*15 ];
  
  for (l = 0; l < 15; l++)
    {
      f_neq[l] = f[l];
    }
  boundary_id = (net_site_data[ i ] & BOUNDARY_ID_MASK) >> BOUNDARY_ID_SHIFT;
  
  *density = inlet_density[ boundary_id ];
  
  *v_x = *v_y = *v_z = 0.F;
  
  f_neq[0] -= (f_new[ f_id[i*15]+cycle_tag*15 ] = f[0] = (2.0/9.0) * *density);
  
  temp = (1.0/9.0) * *density;
  
  for (l = 1; l < 7; l++)
    f_neq[l] -= (f_new[ f_id[i*15+l]+cycle_tag*15 ] = f[l] = temp);
  
  temp *= (1.0/8.0);
  
  for (l = 7; l < 15; l++)
    f_neq[l] -= (f_new[ f_id[i*15+l]+cycle_tag*15 ] = f[l] = temp);
}


// The same as lbmCollision5 but useful for convergence purposes.
void lbmCollisionConv5 (double omega, int i,
			double *density, double *v_x, double *v_y, double *v_z,
			double f_neq[])
{
  double *f;
  double temp;
  
  int l;
  
  unsigned int boundary_id;
  
  
  f = &f_old[ i*30+cycle_tag*15 ];
  
  for (l = 0; l < 15; l++)
    {
      f_neq[l] = f[l];
    }
  boundary_id = (net_site_data[ i ] & BOUNDARY_ID_MASK) >> BOUNDARY_ID_SHIFT;
  
  *density = outlet_density[ boundary_id ];
  
  *v_x = *v_y = *v_z = 0.F;
  
  f_neq[0] -= (f_new[ f_id[i*15]+cycle_tag*15 ] = f[0] = (2.0/9.0) * *density);
  
  temp = (1.0/9.0) * *density;
  
  for (l = 1; l < 7; l++)
    f_neq[l] -= (f_new[ f_id[i*15+l]+cycle_tag*15 ] = f[l] = temp);
  
  temp *= (1.0/8.0);
  
  for (l = 7; l < 15; l++)
    f_neq[l] -= (f_new[ f_id[i*15+l]+cycle_tag*15 ] = f[l] = temp);
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


// von Mises stress computation given the non-equilibrium distribution functions.
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


// Compute the shear stress, i.e. the magnitude of force per unit area
// tangential to the surface with normal "nor[]".
void lbmStress (double density, double f[], double nor[], double *stress)
{
  int e[] = {
    0, 0, 0,
    1, 0, 0,
   -1, 0, 0,
    0, 1, 0,
    0,-1, 0,
    0, 0, 1,
    0, 0,-1,
    1, 1, 1,
   -1,-1,-1,
    1, 1,-1,
   -1,-1, 1,
    1,-1, 1,
   -1, 1,-1,
    1,-1,-1,
   -1, 1, 1
  };
  double sigma[9];                            // stress tensor;
                                              // sigma_ij is the force
                                              // per unit area in
                                              // direction i on the
                                              // plane with the normal
                                              // in direction j
  double stress_vector[] = {0.0, 0.0, 0.0};   // Force per unit area in
                                              // direction i on the
                                              // plane perpendicular to
                                              // the surface normal
  double square_stress_vector = 0.0;
  double normal_stress = 0.0;                 // Magnitude of force per
                                              // unit area normal to the
                                              // surface
  int i, j, l;
  
  
  double temp = lbm_stress_par * (-sqrt(2.0));
  
  for (i = 0; i < 3; i++)
    for (j = 0; j <= i; j++)
      {
	sigma[i*3+j] = 0.0;
	
	for (l = 0; l < 15; l++)
	  {
	    sigma[i*3+j] += f[l] * (e[l*3+i] * e[l*3+j]);
	  }
	sigma[i*3+j] *= temp;
      }
  for (i = 0; i < 3; i++)
    for (j = 0; j < i; j++)
      sigma[j*3+i] = sigma[i*3+j];
  
  for (i = 0; i < 3; i++)
    {
      for (j = 0; j < 3; j++)
	stress_vector[i] += sigma[i*3+j] * nor[j];
      
      square_stress_vector += stress_vector[i] * stress_vector[i];
      normal_stress        += stress_vector[i] * nor[i];
    }
  // shear_stress^2 + normal_stress^2 = stress_vector^2
  *stress = sqrt(square_stress_vector - normal_stress * normal_stress);
}


// Set up of min/max values at the beginning of each pulsatile cycle.
void lbmInitMinMaxValues (void)
{
  lbm_density_min = +1.0e+30;
  lbm_density_max = +1.0e-30;
  
  lbm_velocity_min = +1.0e+30;
  lbm_velocity_max = +1.0e-30;
  
  lbm_stress_min = +1.0e+30;
  lbm_stress_max = +1.0e-30;
}


// Update the min/max values local to the current subdomain.
void lbmUpdateMinMaxValues (double density, double velocity, double stress)
{
  lbm_density_min = (density < lbm_density_min) ? density : lbm_density_min;
  lbm_density_max = (density > lbm_density_max) ? density : lbm_density_max;
  
  lbm_velocity_min = (velocity < lbm_velocity_min) ? velocity : lbm_velocity_min;
  lbm_velocity_max = (velocity > lbm_velocity_max) ? velocity : lbm_velocity_max;
  
  lbm_stress_min = (stress < lbm_stress_min) ? stress : lbm_stress_min;
  lbm_stress_max = (stress > lbm_stress_max) ? stress : lbm_stress_max;
}


// Fluid site updating for benchmarking purposes.
void lbmUpdateSiteDataBench (double omega, int i, double *density, double *vx,double *vy, double *vz, double *velocity,
			     void lbmCollision (double omega, int i,
						double *density, double *v_x, double *v_y, double *v_z,
						double f_neq[]))
{
  double f_neq[15];
  
  lbmCollision (omega, i, density, vx, vy, vz, f_neq);
}


// Fluid site updating for benchmarking plus computation of flow field values for visualisation purposes.
void lbmUpdateSiteDataBenchPlusVis (double omega, int i, double *density, double *vx,double *vy, double *vz, double *velocity,
				    void lbmCollision (double omega, int i,
						       double *density, double *v_x, double *v_y, double *v_z,
						       double f_neq[]))
{
  double f_neq[15];
  double stress;
  
  
  lbmCollision (omega, i, density, vx, vy, vz, f_neq);
  
  *vx *= (1.0 / *density);
  *vy *= (1.0 / *density);
  *vz *= (1.0 / *density);
  *velocity = sqrt(*vx * *vx + *vy * *vy + *vz * *vz);
  
  if (lbm_stress_type == SHEAR_STRESS)
    {
      if (net_site_nor[ i*3 ] >= 1.0e+30)
	{
	  stress = 1.0e+30;
	}
      else
	{
	  lbmStress (*density, f_neq, &net_site_nor[ i*3 ], &stress);
	}
    }
  else
    {
      lbmStress (f_neq, &stress);
    }
  rtUpdateClusterVoxel (i, *density, *velocity, stress);
}


// Fluid site updating for full-production runs.
void lbmUpdateSiteDataSim (double omega, int i, double *density, double *vx,double *vy, double *vz, double *velocity,
			   void lbmCollision (double omega, int i,
					      double *density, double *v_x, double *v_y, double *v_z,
					      double f_neq[]))
{
  double f_neq[15];
  double stress;
  
  
  lbmCollision (omega, i, density, vx, vy, vz, f_neq);
  
  *vx *= (1.0 / *density);
  *vy *= (1.0 / *density);
  *vz *= (1.0 / *density);
  *velocity = sqrt(*vx * *vx + *vy * *vy + *vz * *vz);
  
  if (lbm_stress_type == SHEAR_STRESS)
    {
      if (net_site_nor[ i*3 ] > 1.0e+30)
	{
	  stress = 0.0;
	}
      else
	{
	  lbmStress (*density, f_neq, &net_site_nor[ i*3 ], &stress);
	}
    }
  else
    {
      lbmStress (f_neq, &stress);
    }
  lbmUpdateMinMaxValues (*density, *velocity, stress);
}


// Fluid site updating for full-production runs plus computation of flow field values for visualisation purposes.
void lbmUpdateSiteDataSimPlusVis (double omega, int i, double *density, double *vx,double *vy, double *vz, double *velocity,
				  void lbmCollision (double omega, int i,
						     double *density, double *v_x, double *v_y, double *v_z,
						     double f_neq[]))
{
  double f_neq[15];
  double stress;
  
  
  lbmCollision (omega, i, density, vx, vy, vz, f_neq);
  
  *vx *= (1.0 / *density);
  *vy *= (1.0 / *density);
  *vz *= (1.0 / *density);
  *velocity = sqrt(*vx * *vx + *vy * *vy + *vz * *vz);
  
  if (lbm_stress_type == SHEAR_STRESS)
    {
      if (net_site_nor[ i*3 ] >= 1.0e+30)
	{
	  lbmUpdateMinMaxValues (*density, *velocity, 0.0);
	  rtUpdateClusterVoxel (i, *density, *velocity, 1.0e+30F);
	}
      else
	{
	  lbmStress (*density, f_neq, &net_site_nor[ i*3 ], &stress);
	  lbmUpdateMinMaxValues (*density, *velocity, stress);
	  rtUpdateClusterVoxel (i, *density, *velocity, stress);
	}
    }
  else
    {
      lbmStress (f_neq, &stress);
      lbmUpdateMinMaxValues (*density, *velocity, stress);
      rtUpdateClusterVoxel (i, *density, *velocity, stress);
    }
}


// Returns the type of collision/streaming update for the fluid site
// with data "site_data".
int lbmCollisionType (unsigned int site_data)
{
  unsigned int boundary_type;
  
  
  if (site_data == FLUID_TYPE)
    {
      return FLUID;
    }
  boundary_type = site_data & SITE_TYPE_MASK;
  
  if (boundary_type == FLUID_TYPE)
    {
      return EDGE;
    }
  if (!(site_data & PRESSURE_EDGE_MASK))
    {
      if (boundary_type == INLET_TYPE)
	{
	  return INLET;
	}
      else
	{
	  return OUTLET;
	}
    }
  else
    {
      if (boundary_type == INLET_TYPE)
	{
	  return INLET | EDGE;
	}
      else
	{
	  return OUTLET | EDGE;
	}
    }
}


// Calculate the BCs for each boundary site type and the
// non-equilibrium distribution functions.
void lbmCalculateBC (double f[], unsigned int site_data, double *density,
		     double *vx, double *vy, double *vz, double f_neq[])
{
  double dummy_density;
  double temp;
  
  int unknowns, l;
  
  unsigned int boundary_type, boundary_config, boundary_id;
  
  
  for (l = 0; l < 15; l++)
    {
      f_neq[ l ] = f[ l ];
    }
  boundary_type = site_data & SITE_TYPE_MASK;
  
  if (boundary_type == FLUID_TYPE)
    {
      *density = 0.;

      for (l = 0; l < 15; l++) *density += f[ l ];
      
      f[0] = (2.0/9.0) * *density;
      
      temp = (1.0/9.0) * *density;
      
      for (l = 1; l < 7; l++) f[l] = temp;
      
      temp *= (1.0/8.0);
      
      for (l = 7; l < 15; l++) f[l] = temp;
      
      *vx = *vy = *vz = 0.F;
    }
  else
    {
      boundary_config = (site_data & BOUNDARY_CONFIG_MASK) >> BOUNDARY_CONFIG_SHIFT;
      boundary_id     = (site_data & BOUNDARY_ID_MASK)     >> BOUNDARY_ID_SHIFT;
      
      unknowns = 0;
      
      for (l = 0; l < 14; l++)
	{
	  if (!(boundary_config & (1U << l))) ++unknowns;
	}
      if (boundary_type == INLET_TYPE)
	{
	  *density = inlet_density[ boundary_id ];
	}
      else
	{
	  *density = outlet_density[ boundary_id ];
	}
      if (unknowns <= 1000000)
	{
	  lbmDensityAndVelocity (f, &dummy_density, vx, vy, vz);
	  lbmFeq (*density, *vx, *vy, *vz, f);
	}
      else
	{
	  f[0] = (2.0/9.0) * *density;
	  
	  temp = (1.0/9.0) * *density;
	  
	  for (l = 1; l < 7; l++) f[l] = temp;
	  
	  temp *= (1.0/8.0);
	  
	  for (l = 7; l < 15; l++) f[l] = temp;
	  
	  *vx = *vy = *vz = 0.F;
	}
    }
  for (l = 0; l < 15; l++)
    {
      f_neq[ l ] -= f[ l ];
    }
}


void lbmUpdateBoundaryDensities (int cycle_id, int time_step, LBM *lbm)
{
  double w = 2.0 * PI / lbm->period;
  
  for (int i = 0; i < lbm->inlets; i++)
    {
      /*
      double coef[]={434.661,-239.217,28.9842,0.810304,5.88148,-37.8293,-32.4343,33.1995,-25.3035};
      double inlet_pressure = coef[0];
      
      for (int l = 1; l <= 4; l++)
	inlet_pressure += coef[l] * cos(l * w * (double)time_step);
      for (int l = 5; l < 9; l++)
	inlet_pressure += coef[l] * sin((l-4) * w * (double)time_step);
      
      inlet_pressure = inlet_pressure / mmHg_TO_PASCAL + 90;
      
      //if (cycle_id == 1)
      //	{
      //	  double t = time_step / lbm->period;
      //	  
      //	  inlet_pressure = (1.0 - t) * 90 + t * inlet_pressure;
      //	}
      inlet_density[i] = lbmConvertPressureToLatticeUnits (inlet_pressure, lbm) / Cs2;
      */
      inlet_density[i] = inlet_density_avg[i] + inlet_density_amp[i] * cos(w * (double)time_step + inlet_density_phs[i]);
    }
  for (int i = 0; i < lbm->outlets; i++)
    {
      outlet_density[i] = outlet_density_avg[i] + outlet_density_amp[i] * cos(w * (double)time_step + outlet_density_phs[i]);
    }
}


void lbmInit (char *system_file_name, LBM *lbm, Net *net)
{
  lbm->system_file_name = system_file_name;
  
  if (!check_conv)
    {
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
    }
  else
    {
      lbmInnerCollision[0] = lbmInnerCollisionConv0;
      lbmInnerCollision[1] = lbmCollisionConv1;
      lbmInnerCollision[2] = lbmCollisionConv2;
      lbmInnerCollision[3] = lbmCollisionConv3;
      lbmInnerCollision[4] = lbmCollisionConv4;
      lbmInnerCollision[5] = lbmCollisionConv5;
      
      lbmInterCollision[0] = lbmInterCollisionConv0;
      lbmInterCollision[1] = lbmCollisionConv1;
      lbmInterCollision[2] = lbmCollisionConv2;
      lbmInterCollision[3] = lbmCollisionConv3;
      lbmInterCollision[4] = lbmCollisionConv4;
      lbmInterCollision[5] = lbmCollisionConv5;
    }
  lbmUpdateSiteData[0][0] = lbmUpdateSiteDataSim;
  lbmUpdateSiteData[0][1] = lbmUpdateSiteDataSimPlusVis;
  lbmUpdateSiteData[1][0] = lbmUpdateSiteDataBench;
  lbmUpdateSiteData[1][1] = lbmUpdateSiteDataBenchPlusVis;
  
  lbm_terminate_simulation = 0;
}


void lbmSetInitialConditions (LBM *lbm, Net *net)
{
  double *f_old_p, *f_new_p, f_eq[15];
  double density;
  double temp;
  
  int i, l;
  
  
  density = 0.;
  
  for (i = 0; i < lbm->outlets; i++)
    {
      density += outlet_density_avg[i] - outlet_density_amp[i];
    }
  density /= lbm->outlets;
  
  for (i = 0; i < net->my_sites; i++)
    {
      f_eq[ 0 ] = (2.0/9.0) * density;
      
      temp = (1.0/9.0) * density;
      
      for (l = 1; l < 7; l++) f_eq[ l ] = temp;
      
      temp *= (1.0/8.0);
      
      for (l = 7; l < 15; l++) f_eq[ l ] = temp;
      
      if (!check_conv)
	{
	  f_old_p = &f_old[ i*15 ];
	  f_new_p = &f_new[ i*15 ];
	  
	  for (l = 0; l < 15; l++)
	    {
	      f_new_p[ l ] = f_old_p[ l ] = f_eq[ l ];
	    }
	}
      else
	{
	  for (cycle_tag = 0; cycle_tag < 2; cycle_tag++)
	    {
	      f_old_p = &f_old[ (i*2+cycle_tag)*15 ];
	      f_new_p = &f_new[ (i*2+cycle_tag)*15 ];
	      
	      for (l = 0; l < 15; l++)
		{
		  f_new_p[ l ] = f_old_p[ l ] = f_eq[ l ];
		}
	    }
	}
    }
}


// The entire simulation time step takes place through this function
// when the convergence criterion is not applied. Communications
// automatically handle the streaming stage pertaining to neighbouring
// subdomains.
int lbmCycle (int perform_rt, LBM *lbm, Net *net)
{
  double omega;
  double density, vx, vy, vz, velocity;
  
  int collision_type, offset;
  int i, m;
  
  NeighProc *neigh_proc_p;
  
  
  for (m = 0; m < net->neigh_procs; m++)
    {
      neigh_proc_p = &net->neigh_proc[ m ];
#ifndef NOMPI
      net->err = MPI_Irecv (&f_old[ neigh_proc_p->f_head ],
			    neigh_proc_p->fs, MPI_DOUBLE,
			    neigh_proc_p->id, 10, MPI_COMM_WORLD,
			    &net->req[ 0 ][ m ]);
#endif
    }
  
  omega = lbm->omega;
  
  offset = net->my_inner_sites;
  
  for (collision_type = 0; collision_type < COLLISION_TYPES; collision_type++)
    {
      for (i = offset; i < offset + net->my_inter_collisions[ collision_type ]; i++)
	{
	  (*lbmUpdateSiteData[is_bench][perform_rt]) (omega, i, &density, &vx, &vy, &vz, &velocity,
						      lbmInterCollision[ collision_type ]);
	}
      offset += net->my_inter_collisions[ collision_type ];
    }
  
  for (m = 0; m < net->neigh_procs; m++)
    {
      neigh_proc_p = &net->neigh_proc[ m ];
#ifndef NOMPI
      net->err = MPI_Isend (&f_new[ neigh_proc_p->f_head ],
			    neigh_proc_p->fs, MPI_DOUBLE,
			    neigh_proc_p->id, 10, MPI_COMM_WORLD,
			    &net->req[ 0 ][ net->neigh_procs + m ]);
#endif
    }
  
  offset = 0;
  
  for (collision_type = 0; collision_type < COLLISION_TYPES; collision_type++)
    {
      for (i = offset; i < offset + net->my_inner_collisions[ collision_type ]; i++)
	{
	  (*lbmUpdateSiteData[is_bench][perform_rt]) (omega, i, &density, &vx, &vy, &vz, &velocity,
						      lbmInnerCollision[ collision_type ]);
	}
      offset += net->my_inner_collisions[ collision_type ];
    }
  
  for (m = 0; m < net->neigh_procs; m++)
    {
#ifndef NOMPI
      net->err = MPI_Wait (&net->req[ 0 ][ m ], net->status);
      net->err = MPI_Wait (&net->req[ 0 ][ net->neigh_procs + m ], net->status);
#endif
    }
  
  // Copy the distribution functions received from the neighbouring
  // processors into the destination buffer "f_new".
  for (i = 0; i < net->shared_fs; i++)
    {
      f_new[ f_recv_iv[i] ] = f_old[ net->neigh_proc[0].f_head + i ];
    }
  double *temp = f_old;
  f_old = f_new;
  f_new = temp;
  
  return STABLE;
}


// The entire simulation time step takes place through this function
// when the convergence criterion is applied. Communications
// automatically handle the streaming stage pertaining to neighbouring
// subdomains.
int lbmCycle (int cycle_id, int time_step, int perform_rt, LBM *lbm, Net *net)
{
  double omega;
  double density, vx[2], vy[2], vz[2], velocity[2];
  
  double sum1, sum2;
  double local_data[3];
  double global_data[3];
  
  int is_converged;
  int collision_type;
  int offset;
  int i, m;
  
  NeighProc *neigh_proc_p;
  
  
  for (m = 0; m < net->neigh_procs; m++)
    {
      neigh_proc_p = &net->neigh_proc[ m ];
#ifndef NOMPI
      net->err = MPI_Irecv (&neigh_proc_p->f_to_recv[ 0 ],
			    neigh_proc_p->fs * 2, MPI_DOUBLE,
			    neigh_proc_p->id, 10, MPI_COMM_WORLD,
			    &net->req[ 0 ][ m ]);
#endif
    }
  
  is_converged = 0;
  
  if (cycle_id == 1)
    {
      conv_error = 1.0e+30;
      
      vx[0] = 1.0e+30;
      vy[0] = 1.0e+30;
      vz[0] = 1.0e+30;
    }
  else
    {
      if (time_step == 1) conv_error = 0.0;
    }
  sum1 = 0.0;
  sum2 = 0.0;
  
  omega = lbm->omega;
  
  offset = net->my_inner_sites;
  
  for (collision_type = 0; collision_type < COLLISION_TYPES; collision_type++)
    {
      for (i = offset; i < offset + net->my_inter_collisions[ collision_type ]; i++)
	{
	  if (cycle_id != 1)
	    {
	      cycle_tag = 0;
	      (*lbmUpdateSiteData[1][0]) (omega, i, &density, &vx[0], &vy[0], &vz[0], &velocity[0],
					  lbmInterCollision[ collision_type ]);
	      vx[0] *= (1.0 / density);
	      vy[0] *= (1.0 / density);
	      vz[0] *= (1.0 / density);
	    }
	  cycle_tag = 1;
	  (*lbmUpdateSiteData[0][perform_rt]) (omega, i, &density, &vx[1], &vy[1], &vz[1], &velocity[1],
					       lbmInterCollision[ collision_type ]);
	  sum1 += sqrt((vx[1] - vx[0]) * (vx[1] - vx[0]) +
		       (vy[1] - vy[0]) * (vy[1] - vy[0]) +
		       (vz[1] - vz[0]) * (vz[1] - vz[0]));
	  sum2 += velocity[1];
	}
      offset += net->my_inter_collisions[ collision_type ];
    }
  
  // Copy the distribution functions from the source buffer "f_old"
  // into the ones to be sent to neighbouring processors.
  for (m = 0; m < net->neigh_procs; m++)
    {
      neigh_proc_p = &net->neigh_proc[ m ];
      
      for (i = 0; i < neigh_proc_p->fs; i++)
	{
	  neigh_proc_p->f_to_send[ i*2   ] = f_old[ neigh_proc_p->f_send_id[i]    ];
	  neigh_proc_p->f_to_send[ i*2+1 ] = f_old[ neigh_proc_p->f_send_id[i]+15 ];
	}
    }
  
  for (m = 0; m < net->neigh_procs; m++)
    {
      neigh_proc_p = &net->neigh_proc[ m ];
#ifndef NOMPI
      net->err = MPI_Isend (&neigh_proc_p->f_to_send[ 0 ],
			    neigh_proc_p->fs * 2, MPI_DOUBLE,
			    neigh_proc_p->id, 10, MPI_COMM_WORLD,
			    &net->req[ 0 ][ net->neigh_procs + m ]);
#endif
    }
  
  offset = 0;
  
  for (collision_type = 0; collision_type < COLLISION_TYPES; collision_type++)
    {
      for (i = offset; i < offset + net->my_inner_collisions[ collision_type ]; i++)
	{
	  if (cycle_id != 1)
	    {
	      cycle_tag = 0;
	      (*lbmUpdateSiteData[1][0]) (omega, i, &density, &vx[0], &vy[0], &vz[0], &velocity[0],
					  lbmInnerCollision[ collision_type ]);
	      vx[0] *= (1.0 / density);
	      vy[0] *= (1.0 / density);
	      vz[0] *= (1.0 / density);
	    }
	  cycle_tag = 1;
	  (*lbmUpdateSiteData[0][perform_rt]) (omega, i, &density, &vx[1], &vy[1], &vz[1], &velocity[1],
					       lbmInnerCollision[ collision_type ]);
	  sum1 += sqrt((vx[1] - vx[0]) * (vx[1] - vx[0]) +
		       (vy[1] - vy[0]) * (vy[1] - vy[0]) +
		       (vz[1] - vz[0]) * (vz[1] - vz[0]));
	  sum2 += velocity[1];
	}
      offset += net->my_inner_collisions[ collision_type ];
    }
  
  // stability and convergence error checking
  int is_unstable = 0;
    
  for (i = 0; i < net->my_sites; i++)
    {
      for (m = 0; m < 15; m++)
	{
	  if (f_old[ i*30+15+m ] < 0.)
	    {
	      is_unstable = 1;
	    }
	}
    }
  
  for (m = 0; m < net->neigh_procs; m++)
    {
#ifndef NOMPI
      net->err = MPI_Wait (&net->req[ 0 ][ m ], net->status);
      net->err = MPI_Wait (&net->req[ 0 ][ net->neigh_procs + m ], net->status);
#endif
    }
  
  // Copy the distribution functions received from the neighbouring
  // processors into the destination buffer "f_new".
  for (m = 0; m < net->neigh_procs; m++)
    {
      neigh_proc_p = &net->neigh_proc[ m ];
      
      for (i = 0; i < neigh_proc_p->fs; i++)
	{
	  f_new[ neigh_proc_p->f_recv_iv[i]    ] = neigh_proc_p->f_to_recv[ i*2   ];
	  f_new[ neigh_proc_p->f_recv_iv[i]+15 ] = neigh_proc_p->f_to_recv[ i*2+1 ];
	}
    }
  double *temp = f_old;
  f_old = f_new;
  f_new = temp;
  
  
  // Combine stability and convergence-related values from all the processors.
  local_data[ 0 ] = (double)is_unstable;
  local_data[ 1 ] = sum1;
  local_data[ 2 ] = sum2;
  
#ifndef NOMPI
  net->err = MPI_Allreduce (local_data, global_data, 3,
			    MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
  global_data[0] = local_data[0];
  global_data[1] = local_data[1];
  global_data[2] = local_data[2];
#endif
  
  is_unstable = (global_data[ 0 ] >= 1.0);
  
  if (cycle_id > 1)
    {  
      sum1 = global_data[ 1 ];
      sum2 = global_data[ 2 ];
      
      conv_error += sum1 / sum2;
      
      if (time_step == lbm->period)
	{
	  conv_error /= lbm->period;
	  
	  if (conv_error < TOL)
	    {
	      is_converged = 1;
	    }
	}
    }
  if (is_unstable)
    {
      return UNSTABLE;
    }
  else if (!is_converged)
    {
      return STABLE;
    }
  else
    {
      return STABLE_AND_CONVERGED;
    }
}


void lbmCalculateFlowFieldValues (LBM *lbm)
{
  double *local_data;
  double *global_data;
  
  int i;

  int lMaxInlets = UtilityFunctions::max(6+lbm->inlets,2*lbm->inlets);

  local_data = (double *)malloc(sizeof(double) * lMaxInlets);
  global_data = (double *)malloc(sizeof(double) * lMaxInlets);
  
#ifndef NOMPI
  local_data[0] = lbm_density_min;
  local_data[1] = lbm_velocity_min;
  local_data[2] = lbm_stress_min;
  local_data[3] = lbm_density_max;
  local_data[4] = lbm_velocity_max;
  local_data[5] = lbm_stress_max;
  
  memcpy (&local_data[6], lbm_peak_inlet_velocity, sizeof(double) * lbm->inlets);
  
  MPI_Reduce (&local_data[0], &global_data[0], 3, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce (&local_data[3], &global_data[3], 3+lbm->inlets, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  
  lbm_density_min  = global_data[0];
  lbm_velocity_min = global_data[1];
  lbm_stress_min   = global_data[2];
  lbm_density_max  = global_data[3];
  lbm_velocity_max = global_data[4];
  lbm_stress_max   = global_data[5];
  
  memcpy (lbm_peak_inlet_velocity, &global_data[6], sizeof(double) * lbm->inlets);
  
  for (i = 0; i < lbm->inlets; i++)
    {
      local_data[ i ] = lbm_average_inlet_velocity[ i ];
      local_data[ lbm->inlets+i ] = lbm_inlet_count[ i ];
    }
  MPI_Reduce (local_data, global_data, 2*lbm->inlets, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  
  for (i = 0; i < lbm->inlets; i++)
    {
      lbm_average_inlet_velocity[ i ] = global_data[ i ];
      lbm_inlet_count[ i ] = global_data[ lbm->inlets+i ];
    }
#endif // NOMPI
  
  free(global_data);
  free(local_data);
  
  for (i = 0; i < lbm->inlets; i++)
    {
      lbm_average_inlet_velocity[i] /= lbm_inlet_count[i];
      lbm_average_inlet_velocity[i] = lbmConvertVelocityToPhysicalUnits (lbm_average_inlet_velocity[i], lbm);
      lbm_peak_inlet_velocity[i] = lbmConvertVelocityToPhysicalUnits (lbm_peak_inlet_velocity[i], lbm);
    }
  
  vis_pressure_min = lbmConvertPressureToPhysicalUnits (lbm_density_min * Cs2, lbm);
  vis_pressure_max = lbmConvertPressureToPhysicalUnits (lbm_density_max * Cs2, lbm);
  
  vis_velocity_min = lbmConvertVelocityToPhysicalUnits (lbm_velocity_min, lbm);
  vis_velocity_max = lbmConvertVelocityToPhysicalUnits (lbm_velocity_max, lbm);
  
  vis_stress_min = lbmConvertStressToPhysicalUnits (lbm_stress_min, lbm);
  vis_stress_max = lbmConvertStressToPhysicalUnits (lbm_stress_max, lbm);
  
  vis_period = lbm->period;
  
  vis_inlets = lbm->inlets;
}


int lbmIsUnstable (Net *net)
{
  int is_unstable, stability;
  
  is_unstable = 0;
  
  for (int i = 0; i < net->my_sites; i++)
    {
      for (int l = 0; l < 15; l++)
	{
	  if (f_old[ i*15+l ] < 0.)
	    {
	      is_unstable = 1;
	    }
	}
    }
  
#ifndef NOMPI
  net->err = MPI_Allreduce (&is_unstable, &stability, 1,
			    MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  is_unstable = stability;
#endif
  
  return is_unstable;
}


// Update peak and average inlet velocities local to the current subdomain. 
void lbmUpdateInletVelocities (int time_step, LBM *lbm, Net *net)
{
  double density;
  double vx, vy, vz;
  double velocity;
  
  int offset;
  int inlet_id;
  int i;
  int c1, c2;
  
  
  if (time_step == 1)
    {
      for (i = 0; i < lbm->inlets; i++)
	{
	  lbm_peak_inlet_velocity[ i ] = -1e+30;
	  lbm_average_inlet_velocity[ i ] = 0.;
	  lbm_inlet_count[ i ] = 0;
	}
    }
  if (check_conv)
    {
      c1 = 30;
      c2 = 15;
    }
  else
    {
      c1 = 15;
      c2 = 0;
    }
  offset = net->my_inner_collisions[ 0 ] + net->my_inner_collisions[ 1 ];
  
  for (i = offset; i < offset + net->my_inner_collisions[ 2 ]; i++)
    {
      lbmDensityAndVelocity (&f_old[ i*c1+c2 ], &density, &vx, &vy, &vz);
      
      inlet_id = (net_site_data[i] & BOUNDARY_ID_MASK) >> BOUNDARY_ID_SHIFT;
      
      if (is_inlet_normal_available)
	{
	  vx *= lbm_inlet_normal[ 3*inlet_id+0 ];
	  vy *= lbm_inlet_normal[ 3*inlet_id+1 ];
	  vz *= lbm_inlet_normal[ 3*inlet_id+2 ];
	  
	  velocity = vx * vx + vy * vy + vz * vz;
	  
	  if (velocity > 0.)
	    {
	      velocity = sqrt(velocity) / density;
	    }
	  else
	    {
	      velocity = -sqrt(velocity) / density;
	    }
	}
      else
	{
	  velocity = sqrt(vx * vx + vy * vy + vz * vz) / density;
	}
      lbm_peak_inlet_velocity[ inlet_id ] = fmax(lbm_peak_inlet_velocity[ inlet_id ], velocity);
      lbm_average_inlet_velocity[ inlet_id ] += velocity;
      ++lbm_inlet_count[ inlet_id ];
    }
  offset = net->my_inner_sites + net->my_inter_collisions[ 0 ] + net->my_inter_collisions[ 1 ];
  
  for (i = offset; i < offset + net->my_inter_collisions[ 2 ]; i++)
    {
      lbmDensityAndVelocity (&f_old[ i*c1+c2 ], &density, &vx, &vy, &vz);
      
      inlet_id = (net_site_data[i] & BOUNDARY_ID_MASK) >> BOUNDARY_ID_SHIFT;
      
      if (is_inlet_normal_available)
	{
	  vx *= lbm_inlet_normal[ 3*inlet_id+0 ];
	  vy *= lbm_inlet_normal[ 3*inlet_id+1 ];
	  vz *= lbm_inlet_normal[ 3*inlet_id+2 ];
	  
	  velocity = vx * vx + vy * vy + vz * vz;
	  
	  if (velocity > 0.)
	    {
	      velocity = sqrt(velocity) / density;
	    }
	  else
	    {
	      velocity = -sqrt(velocity) / density;
	    }
	}
      else
	{
	  velocity = sqrt(vx * vx + vy * vy + vz * vz) / density;
	}
      lbm_peak_inlet_velocity[ inlet_id ] = fmax(lbm_peak_inlet_velocity[ inlet_id ], velocity);
      lbm_average_inlet_velocity[ inlet_id ] += velocity;
      ++lbm_inlet_count[ inlet_id ];
    }
}


// In the case of instability, this function restart the simulation
// with twice as many time steps per period and update the parameters
// that depends on this change.
void lbmRestart (LBM *lbm, Net *net)
{
  int i;
  
  for (i = 0; i < lbm->inlets; i++)
    {
      inlet_density_avg[i] = lbmConvertPressureToPhysicalUnits (inlet_density_avg[i] * Cs2, lbm);
      inlet_density_amp[i] = lbmConvertPressureGradToPhysicalUnits (inlet_density_amp[i] * Cs2, lbm);
    }
  for (i = 0; i < lbm->outlets; i++)
    {
      outlet_density_avg[i] = lbmConvertPressureToPhysicalUnits (outlet_density_avg[i] * Cs2, lbm);
      outlet_density_amp[i] = lbmConvertPressureGradToPhysicalUnits (outlet_density_amp[i] * Cs2, lbm);
    }
  lbm->period *= 2;
  
  for (i = 0; i < lbm->inlets; i++)
    {
      inlet_density_avg[i] = lbmConvertPressureToLatticeUnits (inlet_density_avg[i], lbm) / Cs2;
      inlet_density_amp[i] = lbmConvertPressureGradToLatticeUnits (inlet_density_amp[i], lbm) / Cs2;
    }
  for (i = 0; i < lbm->outlets; i++)
    {
      outlet_density_avg[i] = lbmConvertPressureToLatticeUnits (outlet_density_avg[i], lbm) / Cs2;
      outlet_density_amp[i] = lbmConvertPressureGradToLatticeUnits (outlet_density_amp[i], lbm) / Cs2;
    }
  lbm->tau = lbmCalculateTau (lbm);
  
  lbm->viscosity = ((2.0 * lbm->tau - 1.0) / 6.0);
  lbm->omega = -1.0 / lbm->tau;
  
  lbm_stress_par = (1.0 - 1.0 / (2.0 * lbm->tau)) / sqrt(2.0);
  
  lbmSetInitialConditions (lbm, net);
}


void lbmEnd (void)
{
  free(outlet_density_avg);
  free(outlet_density_amp);
  free(outlet_density_phs);
  free(outlet_density);
  
  free(inlet_density_avg);
  free(inlet_density_amp);
  free(inlet_density_phs);
  free(inlet_density);
  
  free(lbm_inlet_count);
  free(lbm_inlet_normal);
  free(lbm_average_inlet_velocity);
  free(lbm_peak_inlet_velocity);
}
