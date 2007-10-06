// In this file, the functions useful for the input/output are reported

#include "config.h"


void lbmReadConfig (LBM *lbm, Net *net)
{
  // this function reads the XDR configuration file but does not store the system
  // and calculate some parameters
  
  FILE *system_config;
  XDR xdr_config;
  
  float *block_density_p = NULL;
  float error;
  float temp;
  
  int blocks, blocks_max;
  int iters;
  
  int i, j, k, ii, jj, kk, m, n;
  int flag;
  int dummy = 0;
  
  unsigned int site_data;
  unsigned int site_i, site_j, site_k;
  unsigned int boundary_id;
  
  struct DensityBlock
  {
    int neigh[26];
    int neighs;
    int boundary_sites;
  };
  
  DensityBlock *density_block = NULL;
  DensityBlock *density_block_p = NULL;
  
  
  system_config = fopen (lbm->system_file_name, "r");
  
  xdrstdio_create (&xdr_config, system_config, XDR_DECODE);
  
  xdr_double (&xdr_config, &lbm->lattice_to_system);
  xdr_int    (&xdr_config, &lbm->blocks_x);
  xdr_int    (&xdr_config, &lbm->blocks_y);
  xdr_int    (&xdr_config, &lbm->blocks_z);
  xdr_int    (&xdr_config, &lbm->block_size);
  
  lbm->sites_x = lbm->blocks_x * lbm->block_size;
  lbm->sites_y = lbm->blocks_y * lbm->block_size;
  lbm->sites_z = lbm->blocks_z * lbm->block_size;
  
  lbm->sites_in_a_block = lbm->block_size * lbm->block_size * lbm->block_size;

  i = lbm->block_size;
  
  lbm->shift = 0;

  while (i > 1)
    {
      i >>= 1;

      ++lbm->shift;
    }

  lbm->blocks = lbm->blocks_x * lbm->blocks_y * lbm->blocks_z;
  
  lbm->fluid_sites_per_block = (short int *)malloc(sizeof(short int) * lbm->blocks);

  lbm->total_fluid_sites = 0;
  
  lbm->site_min_x = 1<<30;
  lbm->site_min_y = 1<<30;
  lbm->site_min_z = 1<<30;
  lbm->site_max_x = 0;
  lbm->site_max_y = 0;
  lbm->site_max_z = 0;
  
  blocks_max = dummy;
  blocks = dummy;
  
  if (!lbm->is_checkpoint)
    {
      blocks_max = 10000;
      
      lbm->block_density = (float *)malloc(sizeof(float) * blocks_max);
      
      lbm->block_map = (int *)malloc(sizeof(int) * lbm->blocks);
  
      density_block = (DensityBlock *)malloc(sizeof(DensityBlock) * blocks_max);
    }
  n = -1;
  
  for (i = 0; i < lbm->blocks_x; i++)
    {
      for (j = 0; j < lbm->blocks_y; j++)
	{
	  for (k = 0; k < lbm->blocks_z; k++)
	    {
	      lbm->fluid_sites_per_block[ ++n ] = 0;
	      
	      xdr_int (&xdr_config, &flag);
	      
	      if (!lbm->is_checkpoint) lbm->block_map[ n ] = -1;
	      
	      if (flag == 0) continue;
	      
	      if (!lbm->is_checkpoint)
		{
		  lbm->block_map[ n ] = blocks;
		  
		  if (blocks == blocks_max)
		    {
		      blocks_max <<= 2;
		      lbm->block_density = (float *)realloc(lbm->block_density,
							    sizeof(float) * blocks_max);
		      density_block = (DensityBlock *)realloc(density_block,
							      sizeof(DensityBlock) * blocks_max);
		    }
		  block_density_p = &lbm->block_density[ blocks ];
		  *block_density_p = 0.F;
		  
		  density_block_p = &density_block[ blocks ];
		  density_block_p->boundary_sites = 0;
		  ++blocks;
		}
	      
	      for (ii = 0; ii < lbm->block_size; ii++)
		{
		  site_i = (i << lbm->shift) + ii;
		  
		  for (jj = 0; jj < lbm->block_size; jj++)
		    {
		      site_j = (j << lbm->shift) + jj;
		      
		      for (kk = 0; kk < lbm->block_size; kk++)
			{
			  site_k = (k << lbm->shift) + kk;
			  
			  xdr_u_int (&xdr_config, &site_data);
			  
			  if ((site_data & SITE_TYPE_MASK) == SOLID_TYPE)
			    {
			      continue;
			    }
			  
			  ++lbm->total_fluid_sites;
			  ++lbm->fluid_sites_per_block[ n ];
			  
			  lbm->site_min_x = min(lbm->site_min_x, site_i);
			  lbm->site_min_y = min(lbm->site_min_y, site_j);
			  lbm->site_min_z = min(lbm->site_min_z, site_k);
			  lbm->site_max_x = max(lbm->site_max_x, site_i);
			  lbm->site_max_y = max(lbm->site_max_y, site_j);
			  lbm->site_max_z = max(lbm->site_max_z, site_k);
			  
			  if (!lbm->is_checkpoint)
			    {
			      if ((site_data & SITE_TYPE_MASK) == INLET_TYPE)
				{
				  boundary_id = (site_data & BOUNDARY_ID_MASK) >> BOUNDARY_ID_SHIFT;
				  
				  ++density_block_p->boundary_sites;
				  *block_density_p += lbm->inlet_density[ boundary_id ];
				}
			      else if ((site_data & SITE_TYPE_MASK) == OUTLET_TYPE)
				{
				  boundary_id = (site_data & BOUNDARY_ID_MASK) >> BOUNDARY_ID_SHIFT;
				  
				  ++density_block_p->boundary_sites;
				  *block_density_p += lbm->outlet_density[ boundary_id ];
				}
			    }
			}
		    }
		}
	    }
	}
    }
  xdr_destroy (&xdr_config);
  fclose (system_config);
  
  if (lbm->is_checkpoint) return;
  
  blocks = 0;
  n = -1;
  
  for (i = 0; i < lbm->blocks_x; i++)
    {
      for (j = 0; j < lbm->blocks_y; j++)
	{
	  for (k = 0; k < lbm->blocks_z; k++)
	    {
	      if (lbm->block_map[ ++n ] == -1) continue;
	      
	      block_density_p = &lbm->block_density[ blocks ];
	      
	      density_block_p = &density_block[ blocks ];
	      density_block_p->neighs = 0;
	      ++blocks;
	      
	      if (density_block_p->boundary_sites > 0)
		{
		  *block_density_p /= (float)density_block_p->boundary_sites;
		  continue;
		}
	      for (ii = max(0, i - 1); ii <= min(i + 1, lbm->blocks_x - 1); ii++)
		{
		  for (jj = max(0, j - 1); jj <= min(j + 1, lbm->blocks_y - 1); jj++)
		    {
		      for (kk = max(0, k - 1); kk <= min(k + 1, lbm->blocks_z - 1); kk++)
			{
			  m = (ii * lbm->blocks_y + jj) * lbm->blocks_z + kk;
			  
			  if (m == n || lbm->block_map[ m ] == -1) continue;
			  
			  density_block_p->neigh[ density_block_p->neighs ] = lbm->block_map[ m ];
			  ++density_block_p->neighs;
			}
		    }
		}
	    }
	}
    }
  error = 1.e+30F;
  iters = 0;
  
  while (error > 1.e-3F)
    {
      ++iters;
      
      error = 0.F;
      
      for (n = 0; n < blocks; n++)
	{
	  density_block_p = &density_block[ n ];
	  
	  if (density_block_p->neighs == 0) continue;
	  
	  block_density_p = &lbm->block_density[ n ];
	  
	  temp = *block_density_p;
	  *block_density_p = 0.F;
	  
	  for (m = 0; m < density_block_p->neighs; m++)
	    {
	      *block_density_p += lbm->block_density[ density_block_p->neigh[m] ];
	    }
	  *block_density_p /= (float)density_block_p->neighs;
	  
	  error = fmaxf(error, fabsf(*block_density_p - temp) /
			fmaxf(1.e-30F, *block_density_p));
	}
    }
  free(density_block);
}


#ifdef STEER
void lbmReadParameters (char *parameters_file_name, LBM *lbm, Net *net, SteerParams *steer)
#else
void lbmReadParameters (char *parameters_file_name, LBM *lbm, Net *net)
#endif
{
  // through this function the processor 0 reads the LB parameters
  // and then communicate them to the other processors
  
  double par_to_send[6+2000];
  
  int n;
  int iters;
  
  
  if (net->id == 0)
    {
      FILE *parameters_file = fopen (parameters_file_name, "r");
      
      
      fscanf (parameters_file, "%i\n", &lbm->is_checkpoint);
      fscanf (parameters_file, "%le\n", &lbm->tau);
      fscanf (parameters_file, "%i\n", &lbm->inlets);
      
      if (lbm->inlet_density == NULL)
	{
	  lbm->inlet_density = (double *)malloc(sizeof(double) * lbm->inlets);
	}
      
      for (n = 0; n < lbm->inlets; n++)
	{
	  fscanf (parameters_file, "%le\n", &lbm->inlet_density[ n ]);
	}
      fscanf (parameters_file, "%i\n", &lbm->outlets);
      
      if (lbm->outlet_density == NULL)
	{
	  lbm->outlet_density = (double *)malloc(sizeof(double) * lbm->outlets);
	}
      for (n = 0; n < lbm->outlets; n++)
	{
	  fscanf (parameters_file, "%le\n", &lbm->outlet_density[ n ]);
	}
      fscanf (parameters_file, "%i\n", &lbm->time_steps_max);
      fscanf (parameters_file, "%le\n", &lbm->tolerance);
      fscanf (parameters_file, "%i\n", &lbm->checkpoint_frequency);
      fscanf (parameters_file, "%i\n", &lbm->convergence_frequency);
      
      fclose (parameters_file);
      
      par_to_send[ 0 ] = 0.1 + (double)lbm->inlets;
      par_to_send[ 1 ] = 0.1 + (double)lbm->outlets;
    }
  net->err = MPI_Bcast (par_to_send, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  if (net->id != 0)
    {
      lbm->inlets  = (int)par_to_send[ 0 ];
      lbm->outlets = (int)par_to_send[ 1 ];
      
      if (lbm->inlet_density == NULL)
	{
	  lbm->inlet_density = (double *)malloc(sizeof(double) * lbm->inlets);
	}
      if (lbm->outlet_density == NULL)
	{
	  lbm->outlet_density = (double *)malloc(sizeof(double) * lbm->outlets);
	}
    }
  else
    {
      par_to_send[ 0 ] = 0.1 + (double)lbm->is_checkpoint;
      par_to_send[ 1 ] = lbm->tau;
      
      for (n = 0; n < lbm->inlets; n++)
	{
	  par_to_send[ 2 + n ] = lbm->inlet_density[ n ];
	}
      for (n = 0; n < lbm->outlets; n++)
	{
	  par_to_send[ 2 + lbm->inlets + n ] = lbm->outlet_density[ n ];
	}
      par_to_send[ 2 + lbm->inlets + lbm->outlets ] = 0.1 + (double)lbm->time_steps_max;
      par_to_send[ 3 + lbm->inlets + lbm->outlets ] = lbm->tolerance;
      par_to_send[ 4 + lbm->inlets + lbm->outlets ] = 0.1 + (double)lbm->checkpoint_frequency;
      par_to_send[ 5 + lbm->inlets + lbm->outlets ] = 0.1 + (double)lbm->convergence_frequency;
    }
  net->err = MPI_Bcast (par_to_send, 6 + lbm->inlets + lbm->outlets, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  if (net->id != 0)
    {
      lbm->is_checkpoint = (int)par_to_send[ 0 ];
      lbm->tau           =      par_to_send[ 1 ];
      
      for (n = 0; n < lbm->inlets; n++)
	{
	  lbm->inlet_density[ n ] = par_to_send[ 2 + n ];
	}
      for (n = 0; n < lbm->outlets; n++)
	{
	  lbm->outlet_density[ n ] = par_to_send[ 2 + lbm->inlets + n ];
	}
      lbm->time_steps_max        = (int)par_to_send[ 2 + lbm->inlets + lbm->outlets ];
      lbm->tolerance             =      par_to_send[ 3 + lbm->inlets + lbm->outlets ];
      lbm->checkpoint_frequency  = (int)par_to_send[ 4 + lbm->inlets + lbm->outlets ];
      lbm->convergence_frequency = (int)par_to_send[ 5 + lbm->inlets + lbm->outlets ];
    }
  lbm->viscosity = ((2.0 * lbm->tau - 1.0) / 6.0);
  
  lbm->omega = -1.0 / lbm->tau;
  lbm->stress_par = (1.0 - 1.0 / (2.0 * lbm->tau)) / sqrt(2.0);
  
#ifdef STEER
  steer->tau            = lbm->tau;
  steer->tolerance      = lbm->tolerance;
  steer->max_time_steps = lbm->time_steps_max;
  steer->conv_freq      = lbm->convergence_frequency;
  steer->check_freq     = lbm->checkpoint_frequency;
#endif
}


#ifdef STEER
void lbmUpdateParameters (LBM *lbm, SteerParams *steer)
{
  lbm->tau = steer->tau;
  lbm->tolerance = steer->tolerance;
  lbm->time_steps_max = steer->max_time_steps;
  lbm->convergence_frequency = steer->conv_freq;
  lbm->checkpoint_frequency  = steer->check_freq;
  
  lbm->viscosity = ((2.0 * lbm->tau - 1.0) / 6.0);
  lbm->omega = -1.0 / lbm->tau;
  lbm->stress_par = (1.0 - 1.0 / (2.0 * lbm->tau)) / sqrt(2.0);
}
#endif


void lbmSetInitialConditionsWithCheckpoint (LBM *lbm, Net *net)
{
  // this functions set the initial distribution functions to the
  // equilibrium ones calculated with zero velocity and unitary
  // density if the checkpoint status is zero and accordingly to the
  // flow fields specified in the checkpoint file if the checkpoint
  // status is one
  
  FILE *system_config;
  XDR  xdr_system_config;
  
  double f_eq[15];
  float density;
  float vx, vy, vz;
  float stress;
  
  double *f_old_p, *f_new_p;
  
  int stability;
  int i, l, m, n;
  int flag;
  
  DataBlock *map_block_p;
  
  
  system_config = fopen (lbm->checkpoint_file_name, "r");
  xdrstdio_create (&xdr_system_config, system_config, XDR_DECODE);
  
  xdr_int    (&xdr_system_config, &stability);
  xdr_double (&xdr_system_config, &lbm->lattice_to_system);
  xdr_int    (&xdr_system_config, &lbm->blocks_x);
  xdr_int    (&xdr_system_config, &lbm->blocks_y);
  xdr_int    (&xdr_system_config, &lbm->blocks_z);
  xdr_int    (&xdr_system_config, &lbm->block_size);
  
  for (n = 0; n < lbm->blocks; n++)
    {
      xdr_int (&xdr_system_config, &flag);
      
      if (flag == 0) continue;
      
      map_block_p = &net->map_block[ n ];
      
      for (m = 0; m < net->sites_in_a_block; m++)
	{
	  xdr_int (&xdr_system_config, &flag);
	  
	  if (flag == 0) continue;
	  
	  xdr_float (&xdr_system_config, &density);
	  xdr_float (&xdr_system_config, &vx);
	  xdr_float (&xdr_system_config, &vy);
	  xdr_float (&xdr_system_config, &vz);
	  xdr_float (&xdr_system_config, &stress);
	  
	  if (net->proc_id[ n ] != net->id) continue;
	  
	  lbmFeq (density, vx, vy, vz, f_eq);
	  
	  i = map_block_p->site_data[ m ];
	  
	  f_old_p = &f_old[ i*15 ];
	  f_new_p = &f_new[ i*15 ];
	  
	  for (l = 0; l < 15; l++)
	    {
	      f_new_p[ l ] = f_old_p[ l ] = f_eq[ l ];
	    }
	  vel[ i ].x = vel[ i ].y = vel[ i ].z = 1.e+30;
	}
    }
  xdr_destroy (&xdr_system_config);
  fclose (system_config);
}


void lbmWriteConfig (int stability, char *output_file_name, int is_checkpoint, LBM *lbm, Net *net)
{
  // this function writes the density, velocity and effective von
  // Mises stress flow field on "output_file_name" which will be the
  // output file or the checkpoint one if "is_checkpoint" is "0" or
  // "1" respectively
  
  FILE *system_config = NULL;
  XDR	xdr_system_config;
  
  double density, density_min, density_max;
  double velocity, velocity_min, velocity_max;
  double stress, stress_min, stress_max;
  double vx, vy, vz;
  double f_eq[15], f_neq[15];
  double *f_old_p;
  
  float macroscopic_par_buffer[ MACROSCOPIC_PARS * 16 * 16 * 16 ];
  
  int l, m, n;
  int flag;
  int my_site_id;
  
  unsigned int site_data;
  
  DataBlock *map_block_p;
  
  
  if (net->id == 0)
    {
      system_config = fopen (output_file_name, "w");
      xdrstdio_create (&xdr_system_config, system_config, XDR_ENCODE);
      xdr_int (&xdr_system_config, &stability);
    }
  
  if (stability == UNSTABLE)
    {
      if (net->id == 0)
	{
	  xdr_destroy(&xdr_system_config);
	  fclose (system_config);
	}
      return;
    }
  
  if (net->id == 0)
    {
      xdr_double (&xdr_system_config, &lbm->lattice_to_system);
      xdr_int    (&xdr_system_config, &lbm->blocks_x);
      xdr_int    (&xdr_system_config, &lbm->blocks_y);
      xdr_int    (&xdr_system_config, &lbm->blocks_z);
      xdr_int    (&xdr_system_config, &lbm->block_size);
    }
  
  density_min = +1.e+30;
  density_max = -1.e+30;
  
  velocity_min = +1.e+30;
  velocity_max = -1.e+30;
  
  stress_min = +1.e+30;
  stress_max = -1.e+30;
  
  for (n = 0; n < lbm->blocks; n++)
    {
      if (net->proc_id[ n ] >= 1 << 14)
	{
	  flag = 0;
	}
      else
	{
	  flag = 1;
	}
      if (net->id == 0)
	{
	  xdr_int (&xdr_system_config, &flag);
	}
      if (flag == 0) continue;
      
      if (net->proc_id[ n ] == net->id)
	{
	  map_block_p = &net->map_block[ n ];
	  
	  for (m = 0; m < net->sites_in_a_block; m++)
	    {
	      my_site_id = map_block_p->site_data[ m ];
	      
	      if (my_site_id & (1U << 31U))
		{
		  macroscopic_par_buffer[m*MACROSCOPIC_PARS+0] = -1.e+30F;
		  continue;
		}
	      
	      f_old_p = &f_old[ my_site_id*15 ];
	      site_data = net->site_data[ my_site_id ];
	      
	      if (site_data == FLUID_TYPE)
		{
		  lbmFeq (f_old_p, &density, &vx, &vy, &vz, f_eq);
		}
	      else
		{
		  lbmCalculateBC (f_old_p, site_data, &density, &vx, &vy, &vz, lbm);
		}
	      for (l = 0; l < 15; l++)
		{
		  f_neq[ l ] = f_old_p[ l ] - f_eq[ l ];
		}
	      macroscopic_par_buffer[m*MACROSCOPIC_PARS+0] = (float)density;
	      macroscopic_par_buffer[m*MACROSCOPIC_PARS+1] = (float)vx;
	      macroscopic_par_buffer[m*MACROSCOPIC_PARS+2] = (float)vy;
	      macroscopic_par_buffer[m*MACROSCOPIC_PARS+3] = (float)vz;
	      macroscopic_par_buffer[m*MACROSCOPIC_PARS+4] =
		(float)(lbm->stress_par * sqrt(lbmStress (f_neq)));
	    }
	}
      if (net->id != 0 && net->proc_id[ n ] == net->id)
	{
	  net->err = MPI_Send (&macroscopic_par_buffer[ 0 ],
			       MACROSCOPIC_PARS * net->sites_in_a_block, MPI_REAL,
			       0, 10, MPI_COMM_WORLD);
	}
      else if (net->id == 0 && net->proc_id[ n ] != net->id)
	{
	  net->err = MPI_Recv (&macroscopic_par_buffer[ 0 ],
			       MACROSCOPIC_PARS * net->sites_in_a_block, MPI_REAL,
			       net->proc_id[ n ], 10, MPI_COMM_WORLD, net->status);
	}
      
      net->err = MPI_Barrier (MPI_COMM_WORLD);
      
      if (net->id != 0) continue;
      
      for (m = 0; m < net->sites_in_a_block; m++)
	{
	  if (macroscopic_par_buffer[ m * MACROSCOPIC_PARS ] < 0.F)
	    {
	      flag = 0;
	    }
	  else
	    {
	      flag = 1;
	    }
	  xdr_int (&xdr_system_config, &flag);
	  
	  if (flag == 0) continue;
	  
	  xdr_float (&xdr_system_config, &macroscopic_par_buffer[m*MACROSCOPIC_PARS+0]);
	  xdr_float (&xdr_system_config, &macroscopic_par_buffer[m*MACROSCOPIC_PARS+1]);
	  xdr_float (&xdr_system_config, &macroscopic_par_buffer[m*MACROSCOPIC_PARS+2]);
	  xdr_float (&xdr_system_config, &macroscopic_par_buffer[m*MACROSCOPIC_PARS+3]);
	  xdr_float (&xdr_system_config, &macroscopic_par_buffer[m*MACROSCOPIC_PARS+4]);
	  
	  density = (double)macroscopic_par_buffer[m*MACROSCOPIC_PARS];
	  density_min = fmin(density_min, density);
	  density_max = fmax(density_max, density);
	  vx = (double)macroscopic_par_buffer[m*MACROSCOPIC_PARS+1];
	  vy = (double)macroscopic_par_buffer[m*MACROSCOPIC_PARS+2];
	  vz = (double)macroscopic_par_buffer[m*MACROSCOPIC_PARS+3];
	  velocity = sqrt(vx * vx + vy * vy + vz * vz);
	  velocity_min = fmin(velocity_min, velocity);
	  velocity_max = fmax(velocity_max, velocity);
	  stress = (double)macroscopic_par_buffer[m*MACROSCOPIC_PARS+4];
	  stress_min = fmin(stress_min, stress);
	  stress_max = fmax(stress_max, stress);
	}
    }
  if (net->id == 0)
    {
      xdr_destroy(&xdr_system_config);
      fclose (system_config);
    }
  lbm->density_min = density_min;
  lbm->density_max = density_max;
  
  lbm->velocity_min = velocity_min;
  lbm->velocity_max = velocity_max;
  
  lbm->stress_min = stress_min;
  lbm->stress_max = stress_max;
}
