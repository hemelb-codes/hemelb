// In this file, the functions useful to calculate the equilibrium distribution
// function, momentums, the effective von Mises stress and the boundary conditions
// are reported

#include "config.h"


void lbmFeq (double f[], double *density, double *v_x, double *v_y, double *v_z, double f_eq[])
{
  double v_xx, v_yy, v_zz;
  double temp1, temp2;
  
  
  *v_x = f[ 1 ] + (f[  7 ] + f[  9 ]) + (f[ 11 ] + f[ 13 ]);
  *v_y = f[ 3 ] + (f[ 12 ] + f[ 14 ]);
  *v_z = f[ 5 ] + f[ 10 ];
  
  *density = f[ 0 ] + (f[ 2 ] + f[ 4 ]) + (f[ 6 ] + f[ 8 ]) + *v_x + *v_y + *v_z;
  
  *v_x -= f[ 2 ] + f[ 8 ] + f[ 10 ] + (f[ 12 ] + f[ 14 ]);
  *v_y += (f[ 7 ] + f[ 9 ]) - (f[ 4 ] + f[ 8 ] + f[ 10 ] + (f[ 11 ] + f[ 13 ]));
  *v_z += f[ 7 ] + f[ 11 ] + f[ 14 ] - ((f[ 6 ] + f[ 8 ]) + f[ 9 ] + f[ 12 ] + f[ 13 ]);

  v_xx = *v_x * *v_x;
  v_yy = *v_y * *v_y;
  v_zz = *v_z * *v_z;
  
  temp1 = (1.0 / 8.0) * *density;
  
  temp2 = v_xx + v_yy + v_zz;
  
  f_eq[ 0 ] = temp1 - (1.0 / 3.0) * temp2;
  
  temp1 -= (1.0 / 6.0) * temp2;
  
  f_eq[  1 ] = temp1 + ((1.0 / 3.0) * *v_x) + 0.5 * v_xx;   // (+1, 0, 0)
  f_eq[  2 ] = temp1 - ((1.0 / 3.0) * *v_x) + 0.5 * v_xx;   // (+1, 0, 0)
  
  f_eq[  3 ] = temp1 + ((1.0 / 3.0) * *v_y) + 0.5 * v_yy;   // (0, +1, 0)
  f_eq[  4 ] = temp1 - ((1.0 / 3.0) * *v_y) + 0.5 * v_yy;   // (0, -1, 0)
  
  f_eq[  5 ] = temp1 + ((1.0 / 3.0) * *v_z) + 0.5 * v_zz;   // (0, 0, +1)
  f_eq[  6 ] = temp1 - ((1.0 / 3.0) * *v_z) + 0.5 * v_zz;   // (0, 0, -1)
  
  temp1 *= (1.0 / 8.0);
  
  temp2 = *v_x + *v_y + *v_z;
  
  f_eq[  7 ] = temp1 + ((1.0 / 24.0) * temp2) + (1.0 / 16.0) * temp2 * temp2;   // (+1, +1, +1)
  f_eq[  8 ] = temp1 - ((1.0 / 24.0) * temp2) + (1.0 / 16.0) * temp2 * temp2;   // (-1, -1, -1)
							     
  temp2 = *v_x + *v_y - *v_z;				     
  							     
  f_eq[  9 ] = temp1 + ((1.0 / 24.0) * temp2) + (1.0 / 16.0) * temp2 * temp2;   // (+1, +1, -1)
  f_eq[ 10 ] = temp1 - ((1.0 / 24.0) * temp2) + (1.0 / 16.0) * temp2 * temp2;   // (-1, -1, +1)
							     
  temp2 = *v_x - *v_y + *v_z;				     
  							     
  f_eq[ 11 ] = temp1 + ((1.0 / 24.0) * temp2) + (1.0 / 16.0) * temp2 * temp2;   // (+1, -1, +1)
  f_eq[ 12 ] = temp1 - ((1.0 / 24.0) * temp2) + (1.0 / 16.0) * temp2 * temp2;   // (-1, +1, -1)
							     
  temp2 = *v_x - *v_y - *v_z;				     
  							     
  f_eq[ 13 ] = temp1 + ((1.0 / 24.0) * temp2) + (1.0 / 16.0) * temp2 * temp2;   // (+1, -1, -1)
  f_eq[ 14 ] = temp1 - ((1.0 / 24.0) * temp2) + (1.0 / 16.0) * temp2 * temp2;   // (-1, +1, +1)
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
  
  f_eq[ 0 ] = temp1 - (1.0 / 3.0) * temp2;
  
  temp1 -= (1.0 / 6.0) * temp2;
  
  f_eq[  1 ] = temp1 + ((1.0 / 3.0) * v_x) + 0.5 * v_xx;   // (+1, 0, 0)
  f_eq[  2 ] = temp1 - ((1.0 / 3.0) * v_x) + 0.5 * v_xx;   // (+1, 0, 0)
  
  f_eq[  3 ] = temp1 + ((1.0 / 3.0) * v_y) + 0.5 * v_yy;   // (0, +1, 0)
  f_eq[  4 ] = temp1 - ((1.0 / 3.0) * v_y) + 0.5 * v_yy;   // (0, -1, 0)
  
  f_eq[  5 ] = temp1 + ((1.0 / 3.0) * v_z) + 0.5 * v_zz;   // (0, 0, +1)
  f_eq[  6 ] = temp1 - ((1.0 / 3.0) * v_z) + 0.5 * v_zz;   // (0, 0, -1)
  
  temp1 *= (1.0 / 8.0);
  
  temp2 = v_x + v_y + v_z;
  
  f_eq[  7 ] = temp1 + ((1.0 / 24.0) * temp2) + (1.0 / 16.0) * temp2 * temp2;   // (+1, +1, +1)
  f_eq[  8 ] = temp1 - ((1.0 / 24.0) * temp2) + (1.0 / 16.0) * temp2 * temp2;   // (-1, -1, -1)
							     
  temp2 = v_x + v_y - v_z;				     
  							     
  f_eq[  9 ] = temp1 + ((1.0 / 24.0) * temp2) + (1.0 / 16.0) * temp2 * temp2;   // (+1, +1, -1)
  f_eq[ 10 ] = temp1 - ((1.0 / 24.0) * temp2) + (1.0 / 16.0) * temp2 * temp2;   // (-1, -1, +1)
							     
  temp2 = v_x - v_y + v_z;				     
  							     
  f_eq[ 11 ] = temp1 + ((1.0 / 24.0) * temp2) + (1.0 / 16.0) * temp2 * temp2;   // (+1, -1, +1)
  f_eq[ 12 ] = temp1 - ((1.0 / 24.0) * temp2) + (1.0 / 16.0) * temp2 * temp2;   // (-1, +1, -1)
							     
  temp2 = v_x - v_y - v_z;				     
  							     
  f_eq[ 13 ] = temp1 + ((1.0 / 24.0) * temp2) + (1.0 / 16.0) * temp2 * temp2;   // (+1, -1, -1)
  f_eq[ 14 ] = temp1 - ((1.0 / 24.0) * temp2) + (1.0 / 16.0) * temp2 * temp2;   // (-1, +1, +1)
}


void lbmDensityAndVelocity (double f[], double *density, double *v_x, double *v_y, double *v_z)
{
  *v_x = f[ 1 ] + (f[  7 ] + f[  9 ]) + (f[ 11 ] + f[ 13 ]);
  *v_y = f[ 3 ] + (f[ 12 ] + f[ 14 ]);
  *v_z = f[ 5 ] + f[ 10 ];
  
  *density = f[ 0 ] + (f[ 2 ] + f[ 4 ]) + (f[ 6 ] + f[ 8 ]) + *v_x + *v_y + *v_z;
  
  *v_x -= (f[ 2 ] + f[ 8 ] + f[ 10 ] + (f[ 12 ] + f[ 14 ]));
  *v_y += (f[ 7 ] + f[ 9 ]) - ((f[ 4 ] + f[ 8 ] + f[ 10 ] + (f[ 11 ] + f[ 13 ])));
  *v_z += f[ 7 ] + f[ 11 ] + f[ 14 ] - (((f[ 6 ] + f[ 8 ]) + f[ 9 ] + f[ 12 ] + f[ 13 ]));
}


double lbmStress (double f[])
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
  
  return a + 6.0 * b;
}


void CalculateBC (double f[], unsigned int site_data,double  *vx, double *vy, double *vz, LBM *lbm)
{
  double density, dummy_density;
  
  int unknowns, i;
  
  unsigned int boundary_type, boundary_config, boundary_id;
  
  
  *vx = *vy = *vz = 0.F;
  
  boundary_type = site_data & SITE_TYPE_MASK;
  
  if (boundary_type == FLUID_TYPE)
    {
      density = 0.;

      for (i = 0; i < 15; i++) density += f[ i ];
      
      density *= (1.0 / 8.0);
      
      for (i = 0; i < 7; i++) f[ i ] = density;
      
      density *= (1.0 / 8.0);
      
      for (i = 7; i < 15; i++) f[ i ] = density;
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
	  density = lbm->inlet_density[ boundary_id ];
	}
      else
	{
	  density = lbm->outlet_density[ boundary_id ];
	}
      if (unknowns <= 5)
	{
	  lbmDensityAndVelocity (f, &dummy_density, vx, vy, vz);
	  lbmFeq (density, *vx, *vy, *vz, f);
	}
      else
	{
	  density *= (1.0 / 8.0);
	  
	  for (i = 0; i < 7; i++) f[ i ] = density;
	  
	  density *= (1.0 / 8.0);
	  
	  for (i = 7; i < 15; i++) f[ i ] = density;
	}
    }
}
