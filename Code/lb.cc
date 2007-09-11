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


void lbmVelocity (double f[], double *v_x, double *v_y, double *v_z)
{
  *v_x = f[ 1 ] + (f[  7 ] + f[  9 ]) + (f[ 11 ] + f[ 13 ]);
  *v_y = f[ 3 ] + (f[ 12 ] + f[ 14 ]);
  *v_z = f[ 5 ] + f[ 10 ];
  
  *v_x -= (f[ 2 ] + f[ 8 ] + f[ 10 ] + (f[ 12 ] + f[ 14 ]));
  *v_y += (f[ 7 ] + f[ 9 ]) - ((f[ 4 ] + f[ 8 ] + f[ 10 ] + (f[ 11 ] + f[ 13 ])));
  *v_z += f[ 7 ] + f[ 11 ] + f[ 14 ] - (((f[ 6 ] + f[ 8 ]) + f[ 9 ] + f[ 12 ] + f[ 13 ]));
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


void lbmCalculateBC (double f[], unsigned int site_data, double *density,
		     double *vx, double *vy, double *vz, LBM *lbm)
{
  double dummy_density;
  
  int unknowns, i;
  
  unsigned int boundary_type, boundary_config, boundary_id;
  
  
  *vx = *vy = *vz = 0.F;
  
  boundary_type = site_data & SITE_TYPE_MASK;
  
  if (boundary_type == FLUID_TYPE)
    {
      *density = 0.;

      for (i = 0; i < 15; i++) *density += f[ i ];
      
      *density *= (1.0 / 8.0);
      
      for (i = 0; i < 7; i++) f[ i ] = *density;
      
      *density *= (1.0 / 8.0);
      
      for (i = 7; i < 15; i++) f[ i ] = *density;
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
	  *density = lbm->inlet_density[ boundary_id ];
	}
      else
	{
	  *density = lbm->outlet_density[ boundary_id ];
	}
      if (unknowns <= 5)
	{
	  lbmDensityAndVelocity (f, &dummy_density, vx, vy, vz);
	  lbmFeq (*density, *vx, *vy, *vz, f);
	}
      else
	{
	  *density *= (1.0 / 8.0);
	  
	  for (i = 0; i < 7; i++) f[ i ] = *density;
	  
	  *density *= (1.0 / 8.0);
	  
	  for (i = 7; i < 15; i++) f[ i ] = *density;
	}
    }
}


void lbmInit (char *system_file_name, char *checkpoint_file_name,
	      LBM *lbm, Net *net)
{
  // basically, this function calls other ones only
  
  lbm->system_file_name     = system_file_name;
  lbm->checkpoint_file_name = checkpoint_file_name;
  
  lbmReadConfig (lbm, net);
  lbmSetInitialConditions (lbm, net);
}


void lbmSetOptimizedInitialConditions (LBM *lbm, Net *net)
{
  double *f_old_p, *f_new_p, f_eq[15];
  double density;
  double density_diff;
  
  int neigh_i, neigh_j, neigh_k;
  int i, j, k;
  int ii, jj, kk;
  int l, m, n;
  int sites_a, sites_b;
  int index_a;
  int site_id, global_site_id;
  int iterations;
  int are_fluid_sites_incrementing;
  int first_outlet_iters;
  
  unsigned int *iter_p;
  
  DataBlock *iter_block_p, *map_block_p;
  
  SiteLocation *site_location_a_p, *site_location_b_p;
  
  
  sites_a = lbm->inlet_sites;
  
  iterations = 0;
  
  are_fluid_sites_incrementing = 1;
  
  while (are_fluid_sites_incrementing)
    {
      ++iterations;
      
      sites_b = 0;
      
      are_fluid_sites_incrementing = 0;
      
      for (index_a = 0; index_a < sites_a; index_a++)
	{
	  site_location_a_p = &lbm->site_location_a[ index_a ];
	  
	  for (l = 1; l < 15; l++)
	    {
	      neigh_i = site_location_a_p->i + e_x[ l ];
	      neigh_j = site_location_a_p->j + e_y[ l ];
	      neigh_k = site_location_a_p->k + e_z[ l ];
	      
	      if (neigh_i < 0 || neigh_i >= lbm->sites_x) continue;
	      if (neigh_j < 0 || neigh_j >= lbm->sites_y) continue;
	      if (neigh_k < 0 || neigh_k >= lbm->sites_z) continue;
	      
	      i = neigh_i >> lbm->shift;
	      j = neigh_j >> lbm->shift;
	      k = neigh_k >> lbm->shift;
	      
	      iter_block_p = &lbm->iter_block[(i * lbm->blocks_y + j) * lbm->blocks_z + k];
	      
	      if (iter_block_p->site_data == NULL) continue;
	      
	      ii = neigh_i - (i << lbm->shift);
	      jj = neigh_j - (j << lbm->shift);
	      kk = neigh_k - (k << lbm->shift);
	      
	      site_id = (((ii << lbm->shift) + jj) << lbm->shift) + kk;
	      
	      iter_p = &iter_block_p->site_data[ site_id ];
	      
	      if (*iter_p == (1U << 31U) || *iter_p != (1U << 15U))
		{
		  continue;
		}
	      *iter_p = iterations;
	      
	      global_site_id == (neigh_i * lbm->sites_y + neigh_j) * lbm->sites_z + neigh_k;
	      
	      if (global_site_id == lbm->first_outlet_site_id)
		{
		  first_outlet_iters == iterations;
		}
	      
	      are_fluid_sites_incrementing = 1;
	      
	      site_location_b_p = &lbm->site_location_b[ sites_b ];
	      site_location_b_p->i = neigh_i;
	      site_location_b_p->j = neigh_j;
	      site_location_b_p->k = neigh_k;
	      ++sites_b;
	    }
	}
      site_location_a_p = lbm->site_location_a;
      lbm->site_location_a = lbm->site_location_b;
      lbm->site_location_b = site_location_a_p;
      
      sites_a = sites_b;
    }
  free(lbm->site_location_b);
  free(lbm->site_location_a);
  
  
  density_diff = (lbm->inlet_density[ 0 ] - lbm->outlet_density[ lbm->first_outlet_ref ]) / (double)first_outlet_iters;
  
  for (n = 0; n < lbm->blocks; n++)
    {
      map_block_p = &net->map_block[ n ];
      
      if (map_block_p->site_data == NULL)
	{
	  continue;
	}
      iter_block_p = &lbm->iter_block[ n ];
      
      for (m = 0; m < lbm->sites_in_a_block; m++)
	{
	  site_id = map_block_p->site_data[ m ];
	  
	  if (site_id & (1U << 31U))
	    {
	      continue;
	    }
	  f_old_p = &f_old[ site_id*15 ];
	  f_new_p = &f_new[ site_id*15 ];
	  
	  density = lbm->inlet_density[ 0 ] - density_diff * (double)iter_block_p->site_data[ m ];

	  lbmFeq (density, 0., 0., 0., f_eq);
	  
	  for (l = 0; l < 15; l++)
	    {
	      f_new_p[ l ] = f_old_p[ l ] = f_eq[ l ];
	    }
	  flow_field[ 3 * site_id + 0 ] = (float)density;
	  flow_field[ 3 * site_id + 1 ] = 0.F;
	  flow_field[ 3 * site_id + 2 ] = 0.F;
	}
    }
  for (i = 0; i < net->my_inner_sites + net->my_inter_sites; i++)
    {
      vel[ i ].x = vel[ i ].y = vel[ i ].z = 1.e+30;
    }
  for (n = 0; n < lbm->blocks; n++)
    {
      if (lbm->iter_block[ n ].site_data != NULL)
	{
	  free (lbm->iter_block[ n ].site_data);
	}
    }
  free (lbm->iter_block);
}


int lbmCycle (int write_checkpoint, int check_convergence, int perform_rt, int *is_converged, LBM *lbm, Net *net)
{
  // the entire simulation time step takes place through this function
  
  double seconds;
  double f_eq[15], f_neq[15];
  double omega;
  double density, stress;
  double vx, vy, vz;
  double sum1, sum2;
  double stability_and_convergence_partial[3];
  double stability_and_convergence_total[3];
  double *f_old_p;
  
  int i, l, m, n;
  int is_unstable;
  int *f_id_p;
  
  unsigned int site_data;
  
  Velocity *vel_p;
  
  NeighProc *neigh_proc_p;
  
  
  is_unstable = 0;
  
  *is_converged = 0;
  
  sum1 = 0.0;
  sum2 = 0.0;
  
  omega = lbm->omega;
  
  seconds = myClock ();
  
  for (i = net->my_inner_sites;
       i < net->my_inner_sites + net->my_inter_sites; i++)
    {
      site_data = net->site_data[ i ];
      f_old_p = &f_old[ i*15 ];
      
      stress = 0.;
      
      if (site_data == FLUID_TYPE)
	{
	  lbmFeq (f_old_p, &density, &vx, &vy, &vz, f_eq);
	  
	  for (l = 0; l < 15; l++)
	    {
	      f_old_p[l] += omega * (f_neq[l] = f_old_p[l] - f_eq[l]);
	    }
	  if (perform_rt) stress = lbm->stress_par * sqrt(lbmStress (f_neq));
	}
      else
	{
	  lbmCalculateBC (f_old_p, site_data, &density, &vx, &vy, &vz, lbm);
	}
      vel_p = &vel[ i ];
      
      sum1 += fabs(vel_p->x - vx) + fabs(vel_p->y - vy) + fabs(vel_p->z - vz);
      sum2 += fabs(vx) + fabs(vy) + fabs(vz);
      
      vel_p->x = vx;
      vel_p->y = vy;
      vel_p->z = vz;
      
      if (perform_rt)
	{
	  flow_field[ 3 * i + 0 ] = (float)density;
	  flow_field[ 3 * i + 1 ] = (float)sqrt(vx * vx + vy * vy + vz * vz);
	  flow_field[ 3 * i + 2 ] = (float)stress;
	}
    }
  net->timing[0] += myClock () - seconds;
  seconds = myClock ();
  
  for (m = 0; m < net->neigh_procs; m++)
    {
      neigh_proc_p = &net->neigh_proc[ m ];
      
      for (n = 0; n < neigh_proc_p->fs; n++)
	{
	  neigh_proc_p->f_to_send[ n ] = f_old[ neigh_proc_p->f_send_id[n] ];
	}
    }
  for (m = 0; m < net->inter_m_neigh_procs; m++)
    {
      neigh_proc_p = &net->inter_m_neigh_proc[ m ];
      
      for (n = 0; n < neigh_proc_p->fs; n++)
	{
	  neigh_proc_p->f_to_send[ n ] = f_old[ neigh_proc_p->f_send_id[n] ];
	}
    }
  net->timing[1] += myClock () - seconds;
  seconds = myClock ();
  
  for (m = 0; m < net->neigh_procs; m++)
    {
      neigh_proc_p = &net->neigh_proc[ m ];
      
      net->err = MPI_Issend (&neigh_proc_p->f_to_send[ 0 ],
			     neigh_proc_p->fs, MPI_DOUBLE,
			     neigh_proc_p->id, 10, MPI_COMM_WORLD,
			     &net->req[ 0 ][ net->id * net->procs + m ]);
      
      net->err = MPI_Irecv (&neigh_proc_p->f_to_recv[ 0 ],
			    neigh_proc_p->fs, MPI_DOUBLE,
			    neigh_proc_p->id, 10, MPI_COMM_WORLD,
			    &net->req[ 0 ][ (net->id + net->procs) * net->procs + m ]);
    }
  for (m = 0; m < net->inter_m_neigh_procs; m++)
    {
      neigh_proc_p = &net->inter_m_neigh_proc[ m ];
      
      net->err = MPI_Issend (&neigh_proc_p->f_to_send[ 0 ],
			     neigh_proc_p->fs, MPI_DOUBLE,
			     neigh_proc_p->id, 10, MPI_COMM_WORLD,
			     &net->req[ 1 ][ net->id * net->procs + m ]);
      
      net->err = MPI_Irecv (&neigh_proc_p->f_to_recv[ 0 ],
			    neigh_proc_p->fs, MPI_DOUBLE,
			    neigh_proc_p->id, 10, MPI_COMM_WORLD,
			    &net->req[ 1 ][ (net->id + net->procs) * net->procs + m ]);
    }

  net->timing[2] += myClock () - seconds;
  seconds = myClock ();
  
  for (i = 0; i < net->my_inner_sites; i++)
    {
      site_data = net->site_data[ i ];
      f_old_p = &f_old[ i*15 ];
      
      stress = 0.;
      
      if (site_data == FLUID_TYPE)
	{
	  lbmFeq (f_old_p, &density, &vx, &vy, &vz, f_eq);
	  
	  for (l = 0; l < 15; l++)
	    {
	      f_old_p[l] += omega * (f_neq[ l ] = f_old_p[l] - f_eq[l]);
	    }
	  if (perform_rt) stress = lbm->stress_par * sqrt(lbmStress (f_neq));
	}
      else
	{
	  lbmCalculateBC (f_old_p, site_data, &density, &vx, &vy, &vz, lbm);
	}
      f_id_p = &f_id[ i*15 ];
      
      for (l = 0; l < 15; l++)
	{
	  if (f_old_p[l] < 0.) is_unstable = 1;
	  
	  f_new[ f_id_p[l] ] = f_old_p[l];
	}
      vel_p = &vel[ i ];
      
      sum1 += fabs(vel_p->x - vx) + fabs(vel_p->y - vy) + fabs(vel_p->z - vz);
      sum2 += fabs(vx) + fabs(vy) + fabs(vz);
      
      vel_p->x = vx;
      vel_p->y = vy;
      vel_p->z = vz;
      
      if (perform_rt)
	{
	  flow_field[ 3 * i + 0 ] = (float)density;
	  flow_field[ 3 * i + 1 ] = (float)sqrt(vx * vx + vy * vy + vz * vz);
	  flow_field[ 3 * i + 2 ] = (float)stress;
	}
    }
  net->timing[3] += myClock () - seconds;
  seconds = myClock ();
  
  for (m = 0; m < net->inter_m_neigh_procs; m++)
    {
      net->err = MPI_Wait (&net->req[ 1 ][ net->id * net->procs + m ], net->status);
      net->err = MPI_Wait (&net->req[ 1 ][ (net->id + net->procs) * net->procs + m ], net->status);
    }
  for (m = 0; m < net->neigh_procs; m++)
    {
      net->err = MPI_Wait (&net->req[ 0 ][ net->id * net->procs + m ], net->status);
      net->err = MPI_Wait (&net->req[ 0 ][ (net->id + net->procs) * net->procs + m ], net->status);
    }
  net->timing[2] += myClock () - seconds;
  seconds = myClock ();
  
  for (m = 0; m < net->inter_m_neigh_procs; m++)
    {
      neigh_proc_p = &net->inter_m_neigh_proc[ m ];
      
      for (n = 0; n < neigh_proc_p->fs; n++)
	{
	  f_new[ neigh_proc_p->f_recv_iv[n] ] = neigh_proc_p->f_to_recv[ n ];
	}
    }
  for (m = 0; m < net->neigh_procs; m++)
    {
      neigh_proc_p = &net->neigh_proc[ m ];
      
      for (n = 0; n < neigh_proc_p->fs; n++)
	{
	  f_new[ neigh_proc_p->f_recv_iv[n] ] = neigh_proc_p->f_to_recv[ n ];
	}
    }
  net->timing[1] += myClock () - seconds;
  seconds = myClock ();
  
  for (i = net->my_inner_sites;
       i < net->my_inner_sites + net->my_inter_sites; i++)
    {
      f_old_p = &f_old[ i*15 ];
      f_id_p = &f_id[ i*15 ];
      
      for (l = 0; l < 15; l++)
	{
	  if (f_old_p[l] < 0.) is_unstable = 1;
	  
	  f_new[ f_id_p[l] ] = f_old_p[l];
	}
    }
  f_old_p = f_old;
  f_old = f_new;
  f_new = f_old_p;
  
  net->timing[4] += myClock () - seconds;
  seconds = myClock ();
  
  if (check_convergence)
    {
      if (net->procs > 1)
	{
	  ++net->convergence_count;
	  
	  stability_and_convergence_partial[ 0 ] = (double)is_unstable;
	  stability_and_convergence_partial[ 1 ] = sum1;
	  stability_and_convergence_partial[ 2 ] = sum2;
	  
	  net->err = MPI_Allreduce (stability_and_convergence_partial,
				    stability_and_convergence_total, 3,
				    MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
	  
	  sum1 = stability_and_convergence_total[ 1 ];
	  sum2 = stability_and_convergence_total[ 2 ];
	  
	  is_unstable = (stability_and_convergence_total[ 0 ] >= 1.);
	}
      
      if (sum1 <= sum2 * lbm->tolerance && sum2 > lbm->tolerance)
	{
	  *is_converged = 1;
	}
      if (net->id == 0)
	{
	  printf (" error: %e ", sum1 / sum2);
	}
    }
  if (write_checkpoint)
    {
      lbmWriteConfig (!is_unstable, lbm->checkpoint_file_name, 1, lbm, net);
    }
  
  net->timing[5] += myClock () - seconds;
  
  if (is_unstable)
    {
      return UNSTABLE;
    }
  else
    {
      return STABLE;
    }
}


void lbmEnd (LBM *lbm)
{
  free(lbm->outlet_density);
  lbm->outlet_density = NULL;
  
  free(lbm->inlet_density);
  lbm->inlet_density = NULL;
  
  free(lbm->fluid_sites_per_block);
  lbm->fluid_sites_per_block = NULL;
}
