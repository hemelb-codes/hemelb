// In this file, the functions useful for the input/output are reported

#include "config.h"

void lbmReadConfig (LBM *lbm, Net *net)
{
  // this function reads the XDR configuration file but does not store the system
  // and calculate some parameters
  
  FILE *system_config;
  XDR xdr_config;
  
  int i, j, k, ii, jj, kk, m, n;
  int flag;
  int is_first_outlet_site_found, first_outlet_ref;
  
  unsigned int site_data;
  unsigned int site_i, site_j, site_k;
  unsigned int boundary_id;
  
  SiteLocation *site_location_a_p;
  
  DataBlock *iter_block_p;
  
  
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
  
  
  lbm->iter_block = (DataBlock *)malloc(sizeof(DataBlock) * lbm->blocks);
  
  for (n = 0; n < lbm->blocks; n++)
    {
      lbm->iter_block[ n ].site_data = NULL;
    }
  lbm->inlet_sites = 0;
  
  lbm->site_location_a = (SiteLocation *)malloc(sizeof(SiteLocation) * 1000000);
  lbm->site_location_b = (SiteLocation *)malloc(sizeof(SiteLocation) * 1000000);
  
  lbm->first_outlet_ref = 0;
  
  for (n = 0; n < lbm->outlets; n++)
    {
      // if (lbm->outlet_density[ n ] < 0.999999 * lbm->inlet_density[ 0 ] ||
      if (lbm->outlet_density[ n ] < 4.999999 * lbm->inlet_density[ 0 ] ||
	  lbm->outlet_density[ n ] > 5.000001 * lbm->inlet_density[ 0 ])
//	  lbm->outlet_density[ n ] > 1.000001 * lbm->inlet_density[ 0 ])
	{
	  lbm->first_outlet_ref = n;
	}
    }
  is_first_outlet_site_found = 0;
  
  n = -1;
  
  for (i = 0; i < lbm->blocks_x; i++)
    {
      for (j = 0; j < lbm->blocks_y; j++)
	{
	  for (k = 0; k < lbm->blocks_z; k++)
	    {
	      lbm->fluid_sites_per_block[ ++n ] = 0;
	      
	      xdr_int (&xdr_config, &flag);
	      
	      if (flag == 0) continue;
	      
	      iter_block_p = &lbm->iter_block[ n ];
	      
	      iter_block_p->site_data =
	      	(unsigned int *)malloc(sizeof(unsigned int) * lbm->sites_in_a_block);
	      
	      m = -1;
	      
	      for (ii = 0; ii < lbm->block_size; ii++)
		{
		  site_i = (i << lbm->shift) + ii;
		  
		  for (jj = 0; jj < lbm->block_size; jj++)
		    {
		      site_j = (j << lbm->shift) + jj;
		      
		      for (kk = 0; kk < lbm->block_size; kk++)
			{
			  site_k = (k << lbm->shift) + kk;
			  
			  ++m;
			  xdr_u_int (&xdr_config, &site_data);
			  
			  if ((site_data & SITE_TYPE_MASK) == SOLID_TYPE)
			    {
			      iter_block_p->site_data[ m ] = 1U << 31U;
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
			  
			  if ((site_data & SITE_TYPE_MASK) == INLET_TYPE)
			    {
			      site_location_a_p = &lbm->site_location_a[ lbm->inlet_sites ];
			      site_location_a_p->i = (i << lbm->shift) + ii;
			      site_location_a_p->j = (j << lbm->shift) + jj;
			      site_location_a_p->k = (k << lbm->shift) + kk;
			      ++lbm->inlet_sites;
			      
			      iter_block_p->site_data[ m ] = 0U;
			    }
			  else
			    {
			      iter_block_p->site_data[ m ] = 1U << 15U;
			    }
			  if ((site_data & SITE_TYPE_MASK) == OUTLET_TYPE &&
			      !is_first_outlet_site_found)
			    {
			      boundary_id = (site_data & BOUNDARY_ID_MASK) >> BOUNDARY_ID_SHIFT;
			      
			      if (boundary_id == lbm->first_outlet_ref)
				{
				  is_first_outlet_site_found = 1;
				  lbm->first_outlet_site_id =
				    (site_i * lbm->sites_y + site_j) * lbm->sites_z + site_k;
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
}

/*
#ifdef STEER
void lbmReadParameters (SteerParams *steer, char *parameters_file_name, LBM *lbm, Net *net)
#else
void lbmReadParameters (char *parameters_file_name, LBM *lbm, Net *net)
#endif
{
  // through this function the processor 0 reads the LB parameters
  // and then communicate them to the other processors
  
  double par_to_send[6+2000];
  
  int n;
  
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
//	int iters;
 //     for (n = 0; n < lbm->outlets; n++)
//	{
//	  fscanf (parameters_file, "%i\n", &iters);
//		lbm->outlet_density[ n ] = 5. - iters * 2.e-4;
//	printf("%le\n", lbm->outlet_density[ n ] );
//	} 
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
      if (lbm->inlet_density == NULL)
	{
	  lbm->inlet_density = (double *)malloc(sizeof(double) * lbm->inlets);
	}
      lbm->outlets = (int)par_to_send[ 1 ];
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
	  par_to_send[ 2 + n ] = 0.1 + (double)lbm->inlet_density[ n ];
	}
      for (n = 0; n < lbm->outlets; n++)
	{
	  par_to_send[ 2 + lbm->inlets + n ] = 0.1 + (double)lbm->outlet_density[ n ];
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
	  lbm->inlet_density[ n ] = (int)par_to_send[ 2 + n ];
	}
      for (n = 0; n < lbm->outlets; n++)
	{
	  lbm->outlet_density[ n ] = (int)par_to_send[ 2 + lbm->inlets + n ];
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
*/

#ifdef STEER
void lbmReadParameters (SteerParams *steer, char *parameters_file_name, LBM *lbm, Net *net)
#else
void lbmReadParameters (char *parameters_file_name, LBM *lbm, Net *net)
#endif
{
  // through this function the processor 0 reads the LB parameters
  // and then communicate them to the other processors
	int n;
  
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

void lbmSetInitialConditions (LBM *lbm, Net *net)
{
  // this functions set the initial distribution functions to the
  // equilibrium ones calculated with zero velocity and unitary
  // density if the checkpoint status is zero and accordingly to the
  // flow fields specified in the checkpoint file if the checkpoint
  // status is one
  
  FILE *system_config;
  XDR  xdr_system_config;
  
  double f_eq[15];
  float pressure, density;
  float vx, vy, vz;
  float stress;
  
  double *f_old_p, *f_new_p;
  
  int stability;
  int i, l, m, n;
  int flag;
  
  DataBlock *map_block_p;
  
  
  if (!lbm->is_checkpoint) return;
  
  if (net->id == 0)
    {
      printf("Opening checkpoint file to read: %s\n", lbm->checkpoint_file_name);
      fflush (stdout);
    }
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
	  
	  xdr_float (&xdr_system_config, &pressure);
	  xdr_float (&xdr_system_config, &vx);
	  xdr_float (&xdr_system_config, &vy);
	  xdr_float (&xdr_system_config, &vz);
	  xdr_float (&xdr_system_config, &stress);
	  
	  if (net->proc_id[ n ] != net->id) continue;
	  
	  density = pressure / Cs2;
	  
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
  // this function writes the pressure, velocity and effective von
  // Mises stress flow field on "output_file_name" which will be the
  // output file or the checkpoint one if "is_checkpoint" is "0" or
  // "1" respectively
  
  FILE *system_config = NULL;
  XDR	xdr_system_config;
  
  double seconds;
  double pressure, pressure_min, pressure_max;
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
  
  
  seconds = myClock ();
  
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
  
  pressure_min = +1.e+30;
  pressure_max = -1.e+30;
  
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
		  lbmFeq (f_old_p, &pressure, &vx, &vy, &vz, f_eq);
		}
	      else
		{
		  lbmCalculateBC (f_old_p, site_data, &pressure, &vx, &vy, &vz, lbm);
		}
	      pressure *= Cs2;
	      
	      for (l = 0; l < 15; l++)
		{
		  f_neq[ l ] = f_old_p[ l ] - f_eq[ l ];
		}
	      macroscopic_par_buffer[m*MACROSCOPIC_PARS+0] = (float)pressure;
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
	  
	  pressure = (double)macroscopic_par_buffer[m*MACROSCOPIC_PARS];
	  pressure_min = fmin(pressure_min, pressure);
	  pressure_max = fmax(pressure_max, pressure);
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
      
      if (!is_checkpoint)
	{
	  printf ("pressure_min, max: %le, %le\n", pressure_min, pressure_max);
	  printf ("velocity_min, max: %le, %le\n", velocity_min, velocity_max);
	  printf ("stress_min, max: %le, %le\n", stress_min, stress_max);
	  fflush (stdout);
	}
    }
  net->timing[8] = myClock () - seconds;
}
