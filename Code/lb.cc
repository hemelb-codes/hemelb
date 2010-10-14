// In this file, the functions useful to calculate the equilibrium distribution
// function, momentums, the effective von Mises stress and the boundary conditions
// are reported

#include <math.h>

#include "lb.h"
#include "utilityFunctions.h"
#include "vis/RayTracer.h"

unsigned int getBoundaryConfig(Net* net, int i)
{
  return (net->net_site_data[ i ] & BOUNDARY_CONFIG_MASK) >> BOUNDARY_CONFIG_SHIFT;  
}

void (*lbmPostTimeStep)(double omega, int i, double *density, double *v_x, double *v_y,
  double *v_z, double f_neq[], Net* net) = NULL;
void (*lbmUpdateSiteData[2])(double omega, int i, double *density, double *vx,
  double *vy, double *vz, double *velocity, Net* net,
  hemelb::lb::collisions::Collision* iCollision);

double LBM::lbmConvertPressureToLatticeUnits (double pressure)
{
  return Cs2 + (pressure - REFERENCE_PRESSURE) * mmHg_TO_PASCAL *
    (PULSATILE_PERIOD / (period * voxel_size)) *
    (PULSATILE_PERIOD / (period * voxel_size)) / BLOOD_DENSITY;
}


double LBM::lbmConvertPressureToPhysicalUnits (double pressure)
{
  return REFERENCE_PRESSURE + ((pressure / Cs2 - 1.0) * Cs2) * BLOOD_DENSITY *
    ((period * voxel_size) / PULSATILE_PERIOD) *
    ((period * voxel_size) / PULSATILE_PERIOD) / mmHg_TO_PASCAL;
}


double LBM::lbmConvertPressureGradToLatticeUnits (double pressure_grad)
{
  return pressure_grad * mmHg_TO_PASCAL *
    (PULSATILE_PERIOD / (period * voxel_size)) *
    (PULSATILE_PERIOD / (period * voxel_size)) / BLOOD_DENSITY;
}


double LBM::lbmConvertPressureGradToPhysicalUnits (double pressure_grad)
{
  return pressure_grad * BLOOD_DENSITY *
    ((period * voxel_size) / PULSATILE_PERIOD) *
    ((period * voxel_size) / PULSATILE_PERIOD) / mmHg_TO_PASCAL;
}


double LBM::lbmConvertVelocityToLatticeUnits (double velocity)
{
  return velocity * (((tau-0.5)/3.0) * voxel_size) / (BLOOD_VISCOSITY / BLOOD_DENSITY);
}


double LBM::lbmConvertVelocityToPhysicalUnits (double velocity)
{
  // convert velocity from lattice units to physical units (m/s)
  return velocity * (BLOOD_VISCOSITY / BLOOD_DENSITY) / (((tau-0.5)/3.0) * voxel_size);
}


double LBM::lbmConvertStressToLatticeUnits (double stress)
{
  return stress * (BLOOD_DENSITY / (BLOOD_VISCOSITY * BLOOD_VISCOSITY)) *
    (((tau-0.5)/3.0) * voxel_size) *
    (((tau-0.5)/3.0) * voxel_size);
}


double LBM::lbmConvertStressToPhysicalUnits (double stress)
{
  // convert stress from lattice units to physical units (Pa)
  return stress * BLOOD_VISCOSITY * BLOOD_VISCOSITY /
    (BLOOD_DENSITY * (((tau-0.5)/3.0) * voxel_size) * (((tau-0.5)/3.0) * voxel_size));
}


void LBM::RecalculateTauViscosityOmega ()
{
  tau = 0.5 + (PULSATILE_PERIOD * BLOOD_VISCOSITY / BLOOD_DENSITY) /
    (Cs2 * period * voxel_size * voxel_size);

  viscosity = ((2.0 * tau - 1.0) / 6.0);
  omega = -1.0 / tau;
  lbm_stress_par = (1.0 - 1.0 / (2.0 * tau)) / sqrt(2.0);
}

//TODO: A lot of this commented-out code doesn't exist anywhere else. Or at
// least some of it doesn't. We should move it all into the class hierarchy
// for collision cases and test it all before removing it here. This needs
// to also be done for the convergence case (which in the new system will simply
// mean calling a certain function with a different variables (like &f_new[whatever + cycle_tag * 15] rather than &f_new[whatever])

/*

// Implementation of simple bounce-back
void simpleBounceBack (double omega, int i,
		    double *density, double *v_x, double *v_y, double *v_z,
		    double f_neq[], Net* net)
{
  double *f = &f_old[ i*15 ];

  lbmDensityAndVelocity (f, density, v_x, v_y, v_z);

  double density_1 = 1. / *density;
  double v_xx = *v_x * *v_x;
  double v_yy = *v_y * *v_y;
  double v_zz = *v_z * *v_z;

  // The actual bounce-back lines, including streaming and collision. Basically swap the non-equilibrium components of f in each of the opposing pairs of directions.
  unsigned int boundary_config = getBoundaryConfig(net, i); 

  int lStreamedIndex[15];

  lStreamedIndex[0] =  f_id[i*15+0];

  for (int l = 1; l < 15; l++)
  {
    if (0 == (boundary_config & (1U << (l-1))))
    {
      // No boundary
      lStreamedIndex[l] = f_id[i*15+l];
    }
    else
    {
      // Yes boundary
      // if direction goes into boundary, do bouce-back. I.e. the f_neq at the current point becomes the f_new in the opposing direction.
      lStreamedIndex[l] = i*15 + (((l % 2) == 0) ? l-1 : l+1);
    }
  }

  f_new[ lStreamedIndex[0] ] = f[0] + omega * (f_neq[0] = f[0] - ((2.0/9.0) * *density - (1.0/3.0) * ((v_xx + v_yy + v_zz) * density_1)));
  
  double temp1 = (1.0/9.0) * *density - (1.0/6.0) * ((v_xx + v_yy + v_zz) * density_1);
  
  f_new[ lStreamedIndex[1] ] = f[1] + omega * (f_neq[1] = f[1] - ((temp1 + (0.5 * density_1) * v_xx) + ((1.0/3.0) * *v_x)));   // (+1, 0, 0)
  f_new[ lStreamedIndex[2] ] = f[2] + omega * (f_neq[2] = f[2] - ((temp1 + (0.5 * density_1) * v_xx) - ((1.0/3.0) * *v_x)));   // (+1, 0, 0)
    
  f_new[ lStreamedIndex[3] ] = f[3] + omega * (f_neq[3] = f[3] - ((temp1 + (0.5 * density_1) * v_yy) + ((1.0/3.0) * *v_y)));   // (0, +1, 0)
  f_new[ lStreamedIndex[4] ] = f[4] + omega * (f_neq[4] = f[4] - ((temp1 + (0.5 * density_1) * v_yy) - ((1.0/3.0) * *v_y)));   // (0, -1, 0)
  
  f_new[ lStreamedIndex[5] ] = f[5] + omega * (f_neq[5] = f[5] - ((temp1 + (0.5 * density_1) * v_zz) + ((1.0/3.0) * *v_z)));   // (0, 0, +1)
  f_new[ lStreamedIndex[6] ] = f[6] + omega * (f_neq[6] = f[6] - ((temp1 + (0.5 * density_1) * v_zz) - ((1.0/3.0) * *v_z)));   // (0, 0, -1)
  
  temp1 *= (1.0/8.0);
  
  double temp2 = (*v_x + *v_y) + *v_z;
  
  f_new[ lStreamedIndex[7] ] = f[7] + omega * (f_neq[7] = f[7] - ((temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) + ((1.0/24.0) * temp2)));   // (+1, +1, +1)
  f_new[ lStreamedIndex[8] ] = f[8] + omega * (f_neq[8] = f[8] - ((temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) - ((1.0/24.0) * temp2)));   // (-1, -1, -1)
  
  temp2 = (*v_x + *v_y) - *v_z;
  
  f_new[ lStreamedIndex[9] ] = f[ 9] + omega * (f_neq[ 9] = f[ 9] - ((temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) + ((1.0/24.0) * temp2)));   // (+1, +1, -1)
  f_new[ lStreamedIndex[10] ] = f[10] + omega * (f_neq[10] = f[10] - ((temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) - ((1.0/24.0) * temp2)));   // (-1, -1, +1)
  
  temp2 = (*v_x - *v_y) + *v_z;
  
  f_new[ lStreamedIndex[11] ] = f[11] + omega * (f_neq[11] = f[11] - ((temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) + ((1.0/24.0) * temp2)));   // (+1, -1, +1)
  f_new[ lStreamedIndex[12] ] = f[12] + omega * (f_neq[12] = f[12] - ((temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) - ((1.0/24.0) * temp2)));   // (-1, +1, -1)
  
  temp2 = (*v_x - *v_y) - *v_z;	 
  
  f_new[ lStreamedIndex[13] ] = f[13] + omega * (f_neq[13] = f[13] - ((temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) + ((1.0/24.0) * temp2)));   // (+1, -1, -1)
  f_new[ lStreamedIndex[14] ] = f[14] + omega * (f_neq[14] = f[14] - ((temp1 + ((1.0/16.0) * density_1) * temp2 * temp2) - ((1.0/24.0) * temp2)));   // (-1, +1, +1)
}

// Implementation of BC3 from Chopard 2008
void regularised (double omega, int i,
		    double *density, double *v_x, double *v_y, double *v_z,
		    double f_neq[], Net* net)
{
  double *f = &f_old[ i*15 ];

  // First calculate the density and macro-velocity
  // TEMPORARILY STORE f_eq IN f_neq BUT THE FUNCTION RETURNS f_eq. THIS IS SORTED 
  // OUT IN THE SUBSEQUENT FOR LOOP.
  lbmFeq (f, density, v_x, v_y, v_z, f_neq);

  // To evaluate PI, first let unknown particle populations take value given by bounce-back of off-equilibrium parts 
  // (fi = fiEq + fopp(i) - fopp(i)Eq)
  unsigned int boundary_config = getBoundaryConfig(net, i);
 
  double fTemp[15];

  for(int l =0; l<15; l++)
    fTemp[l] = f[l];

  if (0 != (boundary_config & (1U << (0))))
  {
    fTemp[1] = f[2] + (2.0/3.0) * *v_x;
  }
  if (0 != (boundary_config & (1U << (1))))
  {
    fTemp[2] = f[1] - (2.0/3.0) * *v_x;
  }
  if (0 != (boundary_config & (1U << (2))))
  {
    fTemp[3] = f[4] + (2.0/3.0) * *v_y;
  }
  if (0 != (boundary_config & (1U << (3))))
  {
    fTemp[4] = f[3] - (2.0/3.0) * *v_y;
  }
  if (0 != (boundary_config & (1U << (4))))
  {
    fTemp[5] = f[6] + (2.0/3.0) * *v_z;
  }
  if (0 != (boundary_config & (1U << (5))))
  {
    fTemp[6] = f[5] - (2.0/3.0) * *v_z;
  }
  if (0 != (boundary_config & (1U << (6))))
  {
    fTemp[7] = f[8] + (2.0/24.0) * ((*v_x + *v_y) + *v_z);
  }
  if (0 != (boundary_config & (1U << (7))))
  {
    fTemp[8] = f[7] - (2.0/24.0) * ((*v_x + *v_y) + *v_z);
  }
  if (0 != (boundary_config & (1U << (8))))
  {
    fTemp[9] = f[10] + (2.0/24.0) * ((*v_x + *v_y) - *v_z);
  }
  if (0 != (boundary_config & (1U << (9))))
  {
    fTemp[10] = f[9] - (2.0/24.0) * ((*v_x + *v_y) - *v_z);
  }
  if (0 != (boundary_config & (1U << (10))))
  {
    fTemp[11] = f[12] + (2.0/24.0) * ((*v_x - *v_y) + *v_z);
  }
  if (0 != (boundary_config & (1U << (11))))
  {
    fTemp[12] = f[11] - (2.0/24.0) * ((*v_x - *v_y) + *v_z);
  }
  if (0 != (boundary_config & (1U << (12))))
  {
    fTemp[13] = f[14] + (2.0/24.0) * ((*v_x - *v_y) - *v_z);
  }
  if (0 != (boundary_config & (1U << (13))))
  {
    fTemp[14] = f[13] - (2.0/24.0) * ((*v_x - *v_y) - *v_z);
  }

  // UP TO THIS POINT, F_NEQ ACTUALLY CONTAINS F_EQ. AT THIS
  // STAGE WE REPLACE IT WITH THE ACTUAL NON-EQ VALUE, POST
  // BOUNCING_BACK WHERE NEEDED.
  for(int l = 0; l < 15; l++)
  {
    f_neq[l] = fTemp[l] - f_neq[l];
  }

  double density_1 = 1. / *density;
  double v_xx = *v_x * *v_x;
  double v_yy = *v_y * *v_y;
  double v_zz = *v_z * *v_z;

  // PI = sum_i e_i e_i f_i
  double piMatrix[3][3];

  double diagSum = f_neq[7] + f_neq[8] + f_neq[9] + f_neq[10] + f_neq[11] + f_neq[12] + f_neq[13] + f_neq[14];

  piMatrix[0][0] = f_neq[1] + f_neq[2] + diagSum;
  piMatrix[0][1] = f_neq[7] + f_neq[8] + f_neq[9] + f_neq[10] - (f_neq[11] + f_neq[12] + f_neq[13] + f_neq[14]);
  piMatrix[0][2] = f_neq[7] + f_neq[8]  + f_neq[11] + f_neq[12] - (f_neq[9] + f_neq[10] + f_neq[13] + f_neq[14]);
  piMatrix[1][0] = piMatrix[0][1];
  piMatrix[1][1] = f_neq[3] + f_neq[4] + diagSum;
  piMatrix[1][2] = f_neq[7] + f_neq[8] + f_neq[13] + f_neq[14] - (f_neq[9] + f_neq[10] + f_neq[11] + f_neq[12]);
  piMatrix[2][0] = piMatrix[0][2];
  piMatrix[2][1] = piMatrix[1][2];
  piMatrix[2][2] = f_neq[5] + f_neq[6] + diagSum;

  for(int m=0; m<3; m++)
    for(int n=0; n<3; n++)
      piMatrix[m][n] /= (2.0 * Cs2 * Cs2);

  // Qi = e_i e_i - (speed of sound ^ 2) * Identity
  // Then gi = fiEq + t_i (the 2/9, 1/9, 1/72 stuff) (Qi . PI (inner product)) / 2 * speed of sound^4
  // Or:  gi = fiEq + t_i (the 2/9, 1/9, 1/72 stuff) ((e_i e_i . PI (inner product)) / 2 * speed of sound^4 - specialNumber)
  double specialNumber = (2.0/9.0) * Cs2 * (piMatrix[0][0] + piMatrix[1][1] + piMatrix[2][2]);
  double piMatrixSum = piMatrix[0][0] + piMatrix[0][1] + piMatrix[0][2] + 
                       piMatrix[1][0] + piMatrix[1][1] + piMatrix[1][2] + 
                       piMatrix[2][0] + piMatrix[2][1] + piMatrix[2][2];

  // The gi (here; f) are then collided and streamed
  f_new[ f_id[i*15] ] = ((2.0/9.0) * *density - (1.0/3.0) * ((v_xx + v_yy + v_zz) * density_1)) + (1.0 + omega) * (f_neq[0] = - specialNumber);
  
  double temp1 = (1.0/9.0) * *density - (1.0/6.0) * ((v_xx + v_yy + v_zz) * density_1);
  specialNumber *= 1.0/2.0;

  // Now apply bounce-back to the components that require it, from fTemp
  int lStreamTo[15];
  for(int l = 1; l < 15; l++)
  {
    int oppL = ((l % 2) == 0) ? (l-1) : (l+1);
    if (0 != (boundary_config & (1U << (l-1))))
    {
      lStreamTo[l] = i*15 + oppL;
    }
    else
    {
      lStreamTo[l] = f_id[i*15+l];
    }
  }

  f_new[ lStreamTo[1] ] = temp1 + (0.5 * density_1) * v_xx + (1.0/3.0) * *v_x + (1.0 + omega) * (f_neq[1] = (1.0/9.0) * piMatrix[0][0] - specialNumber);   // (+1, 0, 0)
  f_new[ lStreamTo[2] ] = temp1 + (0.5 * density_1) * v_xx - (1.0/3.0) * *v_x + (1.0 + omega) * (f_neq[2] = (1.0/9.0) * piMatrix[0][0] - specialNumber);   // (+1, 0, 0)
    
  f_new[ lStreamTo[3] ] = temp1 + (0.5 * density_1) * v_yy + (1.0/3.0) * *v_y + (1.0 + omega) * (f_neq[3] = (1.0/9.0) * piMatrix[1][1] - specialNumber);   // (0, +1, 0)
  f_new[ lStreamTo[4] ] = temp1 + (0.5 * density_1) * v_yy - (1.0/3.0) * *v_y + (1.0 + omega) * (f_neq[4] = (1.0/9.0) * piMatrix[1][1] - specialNumber);   // (0, +1, 0)
  
  f_new[ lStreamTo[5] ] = temp1 + (0.5 * density_1) * v_zz + (1.0/3.0) * *v_z + (1.0 + omega) * (f_neq[5] = (1.0/9.0) * piMatrix[2][2] - specialNumber);   // (0, +1, 0)
  f_new[ lStreamTo[6] ] = temp1 + (0.5 * density_1) * v_zz - (1.0/3.0) * *v_z + (1.0 + omega) * (f_neq[6] = (1.0/9.0) * piMatrix[2][2] - specialNumber);   // (0, +1, 0)
  
  temp1 *= (1.0/8.0);
  specialNumber *= (1.0/8.0);

  double temp2 = (*v_x + *v_y) + *v_z;
  
  f_new[ lStreamTo[7] ] = temp1 + (1.0/16.0) * density_1 * temp2 * temp2 + (1.0/24.0) * temp2 + (1.0 + omega) * (f_neq[7] = ((1.0/72.0) * piMatrixSum - specialNumber ));   // (+1, +1, +1)
  f_new[ lStreamTo[8] ] = temp1 + (1.0/16.0) * density_1 * temp2 * temp2 + (-1.0/24.0) * temp2 + (1.0 + omega) * (f_neq[8] = ((1.0/72.0) * piMatrixSum - specialNumber ));   // (-1, -1, -1)
  
  temp2 = (*v_x + *v_y) - *v_z;
  
  f_new[ lStreamTo[9] ] = temp1 + (1.0/16.0) * density_1 * temp2 * temp2 + (1.0/24.0) * temp2 + (1.0 + omega) * (f_neq[9] = ((1.0/72.0) * (piMatrixSum - 4.0 * (piMatrix[0][2] + piMatrix[1][2]))- specialNumber ));   // (+1, +1, -1)
  f_new[ lStreamTo[10] ] = temp1 + (1.0/16.0) * density_1 * temp2 * temp2 + (-1.0/24.0) * temp2 + (1.0 + omega) * (f_neq[10] = ((1.0/72.0) * (piMatrixSum - 4.0 * (piMatrix[0][2] + piMatrix[1][2]))- specialNumber ));   // (-1, -1, +1)
  
  temp2 = (*v_x - *v_y) + *v_z;
  
  f_new[ lStreamTo[11] ] = temp1 + (1.0/16.0) * density_1 * temp2 * temp2 + (1.0/24.0) * temp2 + (1.0 + omega) * (f_neq[11] = ((1.0/72.0) * (piMatrixSum - 4.0 * (piMatrix[0][1] + piMatrix[1][2]))- specialNumber ));   // (+1, -1, +1)
  f_new[ lStreamTo[12] ] = temp1 + (1.0/16.0) * density_1 * temp2 * temp2 + (-1.0/24.0) * temp2 + (1.0 + omega) * (f_neq[12] = ((1.0/72.0) * (piMatrixSum - 4.0 * (piMatrix[0][1] + piMatrix[1][2]))- specialNumber ));   // (-1, +1, -1)
  
  temp2 = (*v_x - *v_y) - *v_z;	 
  
  f_new[ lStreamTo[13] ] = temp1 + (1.0/16.0) * density_1 * temp2 * temp2 + (1.0/24.0) * temp2 + (1.0 + omega) * (f_neq[13] = ((1.0/72.0) * (piMatrixSum - 4.0 * (piMatrix[0][1] + piMatrix[0][2]))- specialNumber ));   // (+1, -1, -1)
  f_new[ lStreamTo[14] ] = temp1 + (1.0/16.0) * density_1 * temp2 * temp2 + (-1.0/24.0) * temp2 + (1.0 + omega) * (f_neq[14] = ((1.0/72.0) * (piMatrixSum - 4.0 * (piMatrix[0][1] + piMatrix[0][2]))- specialNumber ));   // (-1, +1, +1)
}

// Implementation of interpolation of f values based on distance to boundary.
void LBM::fInterpolation (double omega, int i,
		    double *density, double *v_x, double *v_y, double *v_z,
		    double f_neq[], Net* net)
{
  lbmInterCollision0(omega, i, density, v_x, v_y, v_z, f_neq, net);
}

// Implementation of interpolation of f values based on distance to boundary.
void fInterpolationPostStep (double omega, int i,
		    double *density, double *v_x, double *v_y, double *v_z,
		    double f_neq[], Net* net)
{
  // Fill in the unknown distributions

  // Ready to go. net->cut_distances[i*14 + l] is your friend.
  unsigned int boundary_config = getBoundaryConfig(net, i);

  double* cut_dists = &(net->cut_distances[i*14]);

  // Handle odd indices, then evens - it's slightly easier to take the odd
  // and even cases separately.
  for (int l = 1; l < 15; l+= 2)
  {
    // Only do it if there's no boundary in the opposite direction, otherwise, just leave at eqm
    if ((0 != (boundary_config & (1U << (l-1)))))// && (0 == (boundary_config & (1U << (l)))))
    {
      double twoQ = 2.0 * cut_dists[l-1];
      if(twoQ < 1.0)
      {
        f_new[ i*15 + l + 1 ] = f_new[ i*15 + l ] + twoQ * (f_old[ i*15 + l] - f_new[ i*15 + l ]) ;
      }
      else
      {
        f_new[ i*15 + l + 1 ] = f_old[ i*15 + l + 1 ] + (1. / twoQ) * (f_old[ i*15 + l] - f_old[ i*15 + l + 1 ]);
      }
    }
  }

  for (int l = 2; l < 15; l+=2)
  {
    if ((0 != (boundary_config & (1U << (l-1)))))// && (0 == (boundary_config & (1U << (l-2)))))
    {
      double twoQ = 2.0 * cut_dists[l-1];
      if(twoQ < 1.0)
      {
        f_new[ i*15 + l - 1 ] = f_new[ i*15 + l ] + twoQ * (f_old[ i*15 + l] - f_new[ i*15 + l ]) ;
      }
      else
      {
        f_new[ i*15 + l - 1 ] = f_old[ i*15 + l - 1 ] + (1. / twoQ) * (f_old[ i*15 + l] - f_old[ i*15 + l - 1 ]);
      }
    }
  }
}

// Implementation of the Guo, Zheng, Shi boundary condition (2002).
void LBM::gzsBoundary (double omega, int i,
		    double *density, double *v_x, double *v_y, double *v_z,
		    double f_neq[], Net* net)
{
  // First do a normal collision & streaming step, as if we were mid-fluid.
  // NOTE that we use the version that preserves f_old.
  // NOTE that this handily works out the equilibrium density, v_x, v_y and v_z for us
  lbmInnerCollision0(omega, i, density, v_x, v_y, v_z, f_neq, net);

  // Now fill in the un-streamed to distributions (those that point away from boundaries).
  unsigned int boundary_config = getBoundaryConfig(net, i);

  double* cut_dists = &(net->cut_distances[i*14]);

  // Handle odd indices, then evens - it's slightly easier to take the odd
  // and even cases separately.
  for (int l = 1; l < 15; l++)
  {
    if (0 != (boundary_config & (1U << (l-1))))
    {
      int awayFromWallIndex = ((l%2) == 0) ? (l-1) : (l+1);
      double delta = cut_dists[l-1];
      double uWall[3];
      double fNeq;

      // Work out uw1 (noting that ub is 0 until we implement moving walls)
      uWall[0] = (1 - 1./delta) * *v_x;
      uWall[1] = (1 - 1./delta) * *v_y;
      uWall[2] = (1 - 1./delta) * *v_z;
      fNeq = f_neq[awayFromWallIndex];      

      // Interpolate with uw2 if delta < 0.75
      if(delta < 0.75)
      {
      // Only do the extra interpolation if there's gonna be a point there to interpolate from, i.e. there's no boundary
      // in the direction of awayFromWallIndex
        if(0 == (boundary_config & (1U << (awayFromWallIndex - 1))))
        {
          // Need some info about the next node away from the wall in this direction...
          int nextIOut = f_id[i*15 + awayFromWallIndex] / 15;
            double nextNodeDensity, nextNodeV[3], nextNodeFEq[15];
            lbmFeq (&f_old[nextIOut * 15], 
            &nextNodeDensity, 
            &nextNodeV[0], 
            &nextNodeV[1], 
            &nextNodeV[2], 
            &nextNodeFEq[0]);

            for(int a = 0; a < 3; a++)
              uWall[a] = delta * uWall[a] + (1. - delta) * (delta - 1.) * nextNodeV[a] / (1. + delta);

            fNeq = delta * fNeq + (1. - delta) * (f_old[nextIOut*15 + awayFromWallIndex] - nextNodeFEq[awayFromWallIndex]);
        }
        // If there's nothing to extrapolate from we, very lamely, do a 0VE-style operation to fill in the missing velocity.
        else
        {
          for(int a = 0; a < 3; a++)
              uWall[a] = 0.0;//delta * uWall[a];

          fNeq = 0.0;//delta * fNeq;
        }
      }

      // Use a helper function to calculate the actual value of f_eq in the desired direction at the wall node.
      // Note that we assume that the density is the same as at this node
      double fEqTemp[15];
      lbmFeq (*density, uWall[0], uWall[1], uWall[2], fEqTemp);

      // Collide and stream!
      f_new[i*15 + awayFromWallIndex] = fEqTemp[awayFromWallIndex] + (1.0 + omega) * fNeq;
    }
  }
}

// The same as lbmInnerCollision0 but useful for convergence purposes.
void LBM::lbmInnerCollisionConv0 (double omega, int i,
			     double *density, double *v_x, double *v_y, double *v_z,
			     double f_neq[], Net* net)
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
void LBM::lbmInterCollisionConv0 (double omega, int i,
			     double *density, double *v_x, double *v_y, double *v_z,
			     double f_neq[], Net* net)
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
void LBM::lbmCollisionConv1 (double omega, int i,
			double *density, double *v_x, double *v_y, double *v_z,
			double f_neq[], Net* net)
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
void LBM::lbmCollisionConv2 (double omega, int i,
			double *density, double *v_x, double *v_y, double *v_z,
			double f_neq[], Net* net)
{
  double *f;
  double dummy_density;
  
  unsigned int boundary_id, l;
  
  
  f = &f_old[ i*30+cycle_tag*15 ];
  
  for (l = 0; l < 15; l++)
    {
      f_neq[l] = f[l];
    }
  boundary_id = (net->net_site_data[ i ] & BOUNDARY_ID_MASK) >> BOUNDARY_ID_SHIFT;
  
  *density = inlet_density[ boundary_id ];
  
  lbmDensityAndVelocity (f, &dummy_density, v_x, v_y, v_z);
  lbmFeq (*density, *v_x, *v_y, *v_z, f);
  
  for (l = 0; l < 15; l++)
    {
      f_neq[l] -= (f_new[ f_id[i*15+l]+cycle_tag*15 ] = f[l]);
    }
}


// The same as lbmCollision3 but useful for convergence purposes.
void LBM::lbmCollisionConv3 (double omega, int i,
			double *density, double *v_x, double *v_y, double *v_z,
			double f_neq[], Net* net)
{
  double *f;
  double dummy_density;
  
  unsigned int boundary_id, l;
  
  
  f = &f_old[ i*30+cycle_tag*15 ];
  
  for (l = 0; l < 15; l++)
    {
      f_neq[l] = f[l];
    }
  boundary_id = (net->net_site_data[ i ] & BOUNDARY_ID_MASK) >> BOUNDARY_ID_SHIFT;
  
  *density = outlet_density[ boundary_id ];
  
  lbmDensityAndVelocity (f, &dummy_density, v_x, v_y, v_z);
  lbmFeq (*density, *v_x, *v_y, *v_z, f);
  
  for (l = 0; l < 15; l++)
    {
      f_neq[l] -= (f_new[ f_id[i*15+l]+cycle_tag*15 ] = f[l]);
    }
}


// The same as lbmCollision4 but useful for convergence purposes.
void LBM::lbmCollisionConv4 (double omega, int i,
			double *density, double *v_x, double *v_y, double *v_z,
			double f_neq[], Net* net)
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
  boundary_id = (net->net_site_data[ i ] & BOUNDARY_ID_MASK) >> BOUNDARY_ID_SHIFT;
  
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
void LBM::lbmCollisionConv5 (double omega, int i,
			double *density, double *v_x, double *v_y, double *v_z,
			double f_neq[], Net* net)
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
  boundary_id = (net->net_site_data[ i ] & BOUNDARY_ID_MASK) >> BOUNDARY_ID_SHIFT;
  
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
*/

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

// Fluid site updating for full-production runs.
void lbmUpdateSiteDataSim (double omega, int i, double *density, double *vx,double *vy, double *vz, double *velocity, Net* net, 
  hemelb::lb::collisions::Collision* iCollision)
{
  double f_neq[D3Q15::NUMVECTORS];
  double stress;
  
  
  iCollision->DoCollisions(omega, i, density, vx, vy, vz, f_neq, net);
  
  *vx *= (1.0 / *density);
  *vy *= (1.0 / *density);
  *vz *= (1.0 / *density);
  *velocity = sqrt(*vx * *vx + *vy * *vy + *vz * *vz);
  
  if (lbm_stress_type == SHEAR_STRESS)
    {
      if (net->net_site_nor[ i*3 ] > 1.0e+30)
	{
	  stress = 0.0;
	}
      else
	{
	  D3Q15::CalculateShearStress(*density, f_neq, &net->net_site_nor[ i*3 ], &stress, lbm_stress_par);
	}
    }
  else
    {
      D3Q15::CalculateVonMisesStress(f_neq, &stress, lbm_stress_par);
    }
  lbmUpdateMinMaxValues (*density, *velocity, stress);
}


// Fluid site updating for full-production runs plus computation of flow field values for visualisation purposes.
void lbmUpdateSiteDataSimPlusVis (double omega, int i, double *density, double *vx,double *vy, double *vz, double *velocity, Net* net,
  hemelb::lb::collisions::Collision* iCollision)
{
  double f_neq[D3Q15::NUMVECTORS];
  double stress;
  
  
  iCollision->DoCollisions(omega, i, density, vx, vy, vz, f_neq, net);
  
  *vx *= (1.0 / *density);
  *vy *= (1.0 / *density);
  *vz *= (1.0 / *density);
  *velocity = sqrt(*vx * *vx + *vy * *vy + *vz * *vz);
  
  if (lbm_stress_type == SHEAR_STRESS)
    {
      if (net->net_site_nor[ i*3 ] >= 1.0e+30)
	{
	  lbmUpdateMinMaxValues (*density, *velocity, 0.0);
	  hemelb::vis::rtUpdateClusterVoxel (i, *density, *velocity, 1.0e+30F);
	}
      else
	{
	  D3Q15::CalculateShearStress(*density, f_neq, &net->net_site_nor[ i*3 ], &stress, lbm_stress_par);
	  lbmUpdateMinMaxValues (*density, *velocity, stress);
	  hemelb::vis::rtUpdateClusterVoxel (i, *density, *velocity, stress);
	}
    }
  else
    {
      D3Q15::CalculateVonMisesStress(f_neq, &stress, lbm_stress_par);
      lbmUpdateMinMaxValues (*density, *velocity, stress);
      hemelb::vis::rtUpdateClusterVoxel (i, *density, *velocity, stress);
    }
}

// Returns the type of collision/streaming update for the fluid site
// with data "site_data".
unsigned int lbmCollisionType (unsigned int site_data)
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
  
  int l;
  
  unsigned int boundary_type, boundary_id;
  
  for (l = 0; l < D3Q15::NUMVECTORS; l++)
    {
      f_neq[ l ] = f[ l ];
    }
  boundary_type = site_data & SITE_TYPE_MASK;
  
  if (boundary_type == FLUID_TYPE)
    {
      D3Q15::CalculateDensityAndVelocity(f, density, vx, vy, vz);
    }
  else
    {
      boundary_id     = (site_data & BOUNDARY_ID_MASK)     >> BOUNDARY_ID_SHIFT;
      
      if (boundary_type == INLET_TYPE)
	{
	  *density = inlet_density[ boundary_id ];
	}
      else
	{
	  *density = outlet_density[ boundary_id ];
	}

      D3Q15::CalculateDensityAndVelocity(f, &dummy_density, vx, vy, vz);
      D3Q15::CalculateFeq(*density, *vx, *vy, *vz, f);

    }
  for (l = 0; l < D3Q15::NUMVECTORS; l++)
    {
      f_neq[ l ] -= f[ l ];
    }
}


void LBM::lbmUpdateBoundaryDensities (int cycle_id, int time_step)
{
  double w = 2.0 * PI / period;
  
  for (int i = 0; i < inlets; i++)
    {
      inlet_density[i] = inlet_density_avg[i] + inlet_density_amp[i] * cos(w * (double)time_step + inlet_density_phs[i]);
    }
  for (int i = 0; i < outlets; i++)
    {
      outlet_density[i] = outlet_density_avg[i] + outlet_density_amp[i] * cos(w * (double)time_step + outlet_density_phs[i]);
    }
}


void LBM::lbmInit (char *system_file_name_in, char *parameters_file_name, Net *net)
{
  system_file_name = system_file_name_in;

  lbmUpdateSiteData[0] = lbmUpdateSiteDataSim;
  lbmUpdateSiteData[1] = lbmUpdateSiteDataSimPlusVis;

  lbm_terminate_simulation = 0;

  lbmReadConfig (net);

  lbmReadParameters(parameters_file_name, net);

  lbmInitCollisions();
}

void LBM::lbmInitCollisions()
{
  // TODO Note that the convergence checking is not yet implemented in the
  // new boundary condition hierarchy system.
  // It'd be nice to do this with something like
  // MidFluidCollision = new ConvergenceCheckingWrapper(new WhateverMidFluidCollision());

  mMidFluidCollision = new hemelb::lb::collisions::ImplSimpleCollideAndStream();
  mWallCollision = new hemelb::lb::collisions::ImplZeroVelocityEquilibrium();
  mInletCollision
      = new hemelb::lb::collisions::ImplNonZeroVelocityBoundaryDensity(inlet_density);
  mOutletCollision
      = new hemelb::lb::collisions::ImplNonZeroVelocityBoundaryDensity(outlet_density);
  mInletWallCollision
      = new hemelb::lb::collisions::ImplZeroVelocityBoundaryDensity(inlet_density);
  mOutletWallCollision
      = new hemelb::lb::collisions::ImplZeroVelocityBoundaryDensity(outlet_density);

  // TODO: This will eventually be cleverly replaced.
  // Probably with something like a cycle for each boundary condition
  // doing something like Collision->DoPostTimeStep().
  lbmPostTimeStep = NULL;
}

void LBM::lbmSetInitialConditions (Net *net)
{
  double *f_old_p, *f_new_p, f_eq[D3Q15::NUMVECTORS];
  double density;
  double temp;
  
  int i, l;
  
  
  density = 0.;
  
  for (i = 0; i < outlets; i++)
    {
      density += outlet_density_avg[i] - outlet_density_amp[i];
    }
  density /= outlets;
  
  for (i = 0; i < net->my_sites; i++)
    {
    D3Q15::CalculateFeq(density, 0.0, 0.0, 0.0, f_eq);
      
	  f_old_p = &f_old[ i*15 ];
	  f_new_p = &f_new[ i*15 ];
	  
	  for (l = 0; l < D3Q15::NUMVECTORS; l++)
	    {
	      f_new_p[ l ] = f_old_p[ l ] = f_eq[ l ];
	    }
    }
}

// TODO HACK
hemelb::lb::collisions::Collision* LBM::GetCollision(int i)
    {
  switch(i)
  {
    case 0: return mMidFluidCollision;
    case 1: return mWallCollision;
    case 2: return mInletCollision;
    case 3: return mOutletCollision;
    case 4: return mInletWallCollision;
    case 5: return mOutletWallCollision;
  }
  return NULL;
    }

// The entire simulation time step takes place through this function
// when the convergence criterion is not applied. Communications
// automatically handle the streaming stage pertaining to neighbouring
// subdomains.
int LBM::lbmCycle (int perform_rt, Net *net)
{
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
  
  offset = net->my_inner_sites;
  
  for (collision_type = 0; collision_type < COLLISION_TYPES; collision_type++)
  {
    for (i = offset; i < offset + net->my_inter_collisions[collision_type]; i++)
    {
      (*lbmUpdateSiteData[perform_rt]) (omega, i, &density, &vx, &vy, &vz, &velocity, net,
      GetCollision(collision_type));
    }
    offset += net->my_inter_collisions[collision_type];
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
        (*lbmUpdateSiteData[perform_rt]) (omega, i, &density, &vx, &vy, &vz, &velocity, net,
        GetCollision(collision_type));
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

  // Do any cleanup steps necessary on boundary nodes
  // NOTE that 1 is the collision type referring to boundary nodes.
  // One day, this will be a constant.
  if(lbmPostTimeStep != NULL)
  {
    offset = net->my_inner_sites + net->my_inter_collisions[0];
  
    for (i = offset; i < offset + net->my_inter_collisions[ 1 ]; i++)
    {
      (*lbmUpdateSiteData[perform_rt]) (omega, i, &density, &vx, &vy, &vz, &velocity, net,
      GetCollision(collision_type));
    }
  
    offset = net->my_inner_collisions[ 0 ];
  
    for (i = offset; i < offset + net->my_inner_collisions[ 1 ]; i++)
    {
      (*lbmUpdateSiteData[perform_rt]) (omega, i, &density, &vx, &vy, &vz, &velocity, net,
      GetCollision(collision_type));
    }
  }

  // Swap f_old and f_new ready for the next timestep.
  double *temp = f_old;
  f_old = f_new;
  f_new = temp;
  
  return STABLE;
}

void LBM::lbmCalculateFlowFieldValues ()
{
  double *local_data;
  double *global_data;
  
  int i;

  int lMaxInlets = hemelb::util::max(6+inlets,2*inlets);

  local_data = new double[lMaxInlets];
  global_data = new double[lMaxInlets];
  
#ifndef NOMPI
  local_data[0] = lbm_density_min;
  local_data[1] = lbm_velocity_min;
  local_data[2] = lbm_stress_min;
  local_data[3] = lbm_density_max;
  local_data[4] = lbm_velocity_max;
  local_data[5] = lbm_stress_max;
  
  memcpy (&local_data[6], lbm_peak_inlet_velocity, sizeof(double) * inlets);
  
  MPI_Reduce (&local_data[0], &global_data[0], 3, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce (&local_data[3], &global_data[3], 3+inlets, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  
  lbm_density_min  = global_data[0];
  lbm_velocity_min = global_data[1];
  lbm_stress_min   = global_data[2];
  lbm_density_max  = global_data[3];
  lbm_velocity_max = global_data[4];
  lbm_stress_max   = global_data[5];
  
  memcpy (lbm_peak_inlet_velocity, &global_data[6], sizeof(double) * inlets);
  
  for (i = 0; i < inlets; i++)
    {
      local_data[ i ] = lbm_average_inlet_velocity[ i ];
      local_data[ inlets+i ] = lbm_inlet_count[ i ];
    }
  MPI_Reduce (local_data, global_data, 2*inlets, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  
  for (i = 0; i < inlets; i++)
    {
      lbm_average_inlet_velocity[ i ] = global_data[ i ];
      lbm_inlet_count[ i ] = global_data[ inlets+i ];
    }
#endif // NOMPI
  
  delete [] global_data;
  delete [] local_data;
  
  for (i = 0; i < inlets; i++)
    {
      lbm_average_inlet_velocity[i] /= lbm_inlet_count[i];
      lbm_average_inlet_velocity[i] = lbmConvertVelocityToPhysicalUnits (lbm_average_inlet_velocity[i]);
      lbm_peak_inlet_velocity[i] = lbmConvertVelocityToPhysicalUnits (lbm_peak_inlet_velocity[i]);
    }
  
  lbm_phys_pressure_min = lbmConvertPressureToPhysicalUnits (lbm_density_min * Cs2);
  lbm_phys_pressure_max = lbmConvertPressureToPhysicalUnits (lbm_density_max * Cs2);
  
  lbm_phys_velocity_min = lbmConvertVelocityToPhysicalUnits (lbm_velocity_min);
  lbm_phys_velocity_max = lbmConvertVelocityToPhysicalUnits (lbm_velocity_max);
  
  lbm_phys_stress_min = lbmConvertStressToPhysicalUnits (lbm_stress_min);
  lbm_phys_stress_max = lbmConvertStressToPhysicalUnits (lbm_stress_max);
  
  period = period;
  
  inlets = inlets;
}


int lbmIsUnstable (Net *net)
{
  int is_unstable, stability;
  
  is_unstable = 0;
  
  for (int i = 0; i < net->my_sites; i++)
    {
      for (int l = 0; l < D3Q15::NUMVECTORS; l++)
	{
	  if (f_old[ i*D3Q15::NUMVECTORS+l ] < 0.)
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
void LBM::lbmUpdateInletVelocities (int time_step, Net *net)
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
      for (i = 0; i < inlets; i++)
	{
	  lbm_peak_inlet_velocity[ i ] = -1e+30;
	  lbm_average_inlet_velocity[ i ] = 0.;
	  lbm_inlet_count[ i ] = 0;
	}
    }

      c1 = 15;
      c2 = 0;

      offset = net->my_inner_collisions[ 0 ] + net->my_inner_collisions[ 1 ];
  
  for (i = offset; i < offset + net->my_inner_collisions[ 2 ]; i++)
    {
      D3Q15::CalculateDensityAndVelocity(&f_old[ i*c1+c2 ], &density, &vx, &vy, &vz);
      
      inlet_id = (net->net_site_data[i] & BOUNDARY_ID_MASK) >> BOUNDARY_ID_SHIFT;
      
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
      D3Q15::CalculateDensityAndVelocity(&f_old[ i*c1+c2 ], &density, &vx, &vy, &vz);
      
      inlet_id = (net->net_site_data[i] & BOUNDARY_ID_MASK) >> BOUNDARY_ID_SHIFT;
      
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
void LBM::lbmRestart (Net *net)
{
  int i;
  
  for (i = 0; i < inlets; i++)
    {
      inlet_density_avg[i] = lbmConvertPressureToPhysicalUnits (inlet_density_avg[i] * Cs2);
      inlet_density_amp[i] = lbmConvertPressureGradToPhysicalUnits (inlet_density_amp[i] * Cs2);
    }
  for (i = 0; i < outlets; i++)
    {
      outlet_density_avg[i] = lbmConvertPressureToPhysicalUnits (outlet_density_avg[i] * Cs2);
      outlet_density_amp[i] = lbmConvertPressureGradToPhysicalUnits (outlet_density_amp[i] * Cs2);
    }
  period *= 2;
  
  for (i = 0; i < inlets; i++)
    {
      inlet_density_avg[i] = lbmConvertPressureToLatticeUnits (inlet_density_avg[i]) / Cs2;
      inlet_density_amp[i] = lbmConvertPressureGradToLatticeUnits (inlet_density_amp[i]) / Cs2;
    }
  for (i = 0; i < outlets; i++)
    {
      outlet_density_avg[i] = lbmConvertPressureToLatticeUnits (outlet_density_avg[i]) / Cs2;
      outlet_density_amp[i] = lbmConvertPressureGradToLatticeUnits (outlet_density_amp[i]) / Cs2;
    }

  RecalculateTauViscosityOmega ();
  
  lbmSetInitialConditions (net);
}


LBM::~LBM()
{
  // Delete arrays allocated for the inlets
  delete[] inlet_density;
  delete[] inlet_density_avg;
  delete[] inlet_density_amp;
  delete[] inlet_density_phs;

  // Delete arrays allocated for the outlets
  delete[] outlet_density;
  delete[] outlet_density_avg;
  delete[] outlet_density_amp;
  delete[] outlet_density_phs;
  
  // Delete the collision and stream objects we've been using
  delete mMidFluidCollision;
  delete mWallCollision;
  delete mInletCollision;
  delete mOutletCollision;
  delete mInletWallCollision;
  delete mOutletWallCollision;

  // Delete various other arrays used
  delete[] lbm_inlet_count;
  delete[] lbm_inlet_normal;
  delete[] lbm_average_inlet_velocity;
  delete[] lbm_peak_inlet_velocity;
}
