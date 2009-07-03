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
  
  unsigned int site_i, site_j, site_k;
  
  DataBlock *data_block_p;
  ProcBlock *proc_block_p;
  
  
  system_config = fopen (lbm->system_file_name, "r");
  
  if (system_config == NULL) {
    fprintf(stderr, "unable to open file %s [rank %i], exiting\n", lbm->system_file_name, net->id);
    fflush(0x0);
    exit(0x0);
  } else {
    fprintf(stderr, "done\n");
  }
  fflush(NULL);

  xdrstdio_create (&xdr_config, system_config, XDR_DECODE);
  
  xdr_double (&xdr_config, &lbm->lattice_to_system);
  xdr_int    (&xdr_config, &blocks_x);
  xdr_int    (&xdr_config, &blocks_y);
  xdr_int    (&xdr_config, &blocks_z);
  xdr_int    (&xdr_config, &block_size);
  
  sites_x = blocks_x * block_size;
  sites_y = blocks_y * block_size;
  sites_z = blocks_z * block_size;
  
  sites_in_a_block = block_size * block_size * block_size;

  i = block_size;
  
  shift = 0;

  while (i > 1)
    {
      i >>= 1;

      ++shift;
    }
  
  blocks = blocks_x * blocks_y * blocks_z;
  
  net->data_block = (DataBlock *)malloc(sizeof(DataBlock) * blocks);
  
  net->proc_block = (ProcBlock *)malloc(sizeof(ProcBlock) * blocks);
  
  lbm->total_fluid_sites = 0;
  
  lbm->site_min_x = 1<<30;
  lbm->site_min_y = 1<<30;
  lbm->site_min_z = 1<<30;
  lbm->site_max_x = 0;
  lbm->site_max_y = 0;
  lbm->site_max_z = 0;
  
  net->fr_time = myClock ();
  
  n = -1;
  
  for (i = 0; i < blocks_x; i++)
    {
      for (j = 0; j < blocks_y; j++)
	{
	  for (k = 0; k < blocks_z; k++)
	    {
	      ++n;
	      
	      data_block_p = &net->data_block[ n ];
	      proc_block_p = &net->proc_block[ n ];
	      
	      data_block_p->site_data = NULL;
	      proc_block_p->proc_id   = NULL;
	      
	      xdr_int (&xdr_config, &flag);
	      
	      if (flag == 0) continue;
	      
	      data_block_p->site_data = (unsigned int *)malloc(sizeof(unsigned int) * sites_in_a_block);
	      proc_block_p->proc_id   = (int *)malloc(sizeof(int) * sites_in_a_block);
	      
	      m = -1;
	      
	      for (ii = 0; ii < block_size; ii++)
		{
		  site_i = (i << shift) + ii;
		  
		  for (jj = 0; jj < block_size; jj++)
		    {
		      site_j = (j << shift) + jj;
		      
		      for (kk = 0; kk < block_size; kk++)
			{
			  site_k = (k << shift) + kk;
			  
			  ++m;
			  xdr_u_int (&xdr_config, &data_block_p->site_data[ m ]);
			  
			  if ((data_block_p->site_data[ m ] & SITE_TYPE_MASK) == SOLID_TYPE)
			    {
			      proc_block_p->proc_id[ m ] = 1 << 30;
			      continue;
			    }
			  proc_block_p->proc_id[ m ] = -1;
			  
			  ++lbm->total_fluid_sites;
			  
			  lbm->site_min_x = min(lbm->site_min_x, site_i);
			  lbm->site_min_y = min(lbm->site_min_y, site_j);
			  lbm->site_min_z = min(lbm->site_min_z, site_k);
			  lbm->site_max_x = max(lbm->site_max_x, site_i);
			  lbm->site_max_y = max(lbm->site_max_y, site_j);
			  lbm->site_max_z = max(lbm->site_max_z, site_k);
			}
		    }
		}
	    }
	}
    }
  xdr_destroy (&xdr_config);
  fclose (system_config);
  
  net->fr_time = myClock () - net->fr_time;
}


void lbmReadParameters (char *parameters_file_name, LBM *lbm, Net *net)
{
  // through this function the processor 0 reads the LB parameters
  // and then communicate them to the other processors
  
  double par_to_send[10000];
  
  int n;
  
  
  if (net->id == 0)
    {
      FILE *parameters_file = fopen (parameters_file_name, "r");

      if (parameters_file == NULL) {
        fprintf(stderr, "unable to open file %s, exiting\n", parameters_file_name);
        fflush(NULL);
        exit(0x0);
      } else {
        fprintf(stderr, "done\n");
      }
      fflush(NULL);
      
      fscanf (parameters_file, "%i\n", &lbm->inlets);
      
      inlet_density     = (double *)malloc(sizeof(double) * max(1, lbm->inlets));
      inlet_density_avg = (double *)malloc(sizeof(double) * max(1, lbm->inlets));
      inlet_density_amp = (double *)malloc(sizeof(double) * max(1, lbm->inlets));
      inlet_density_phs = (double *)malloc(sizeof(double) * max(1, lbm->inlets));
      
      for (n = 0; n < lbm->inlets; n++)
	{
	  fscanf (parameters_file, "%le %le %le\n",
		  &inlet_density_avg[n], &inlet_density_amp[n], &inlet_density_phs[n]);
	  
	  inlet_density_avg[n] = lbmConvertPressureToLatticeUnits (inlet_density_avg[n], lbm) / Cs2;
	  inlet_density_amp[n] = lbmConvertPressureGradToLatticeUnits (inlet_density_amp[n], lbm) / Cs2;
	  inlet_density_phs[n] *= DEG_TO_RAD;
	  
	  //if (is_bench)
	  //  {
	  //    inlet_density_avg[n] = 1.0;
	  //    inlet_density_amp[n] = 0.0;
	  //    inlet_density_phs[n] = 0.0;
	  //  }
	}
      fscanf (parameters_file, "%i\n", &lbm->outlets);
      
      outlet_density     = (double *)malloc(sizeof(double) * max(1, lbm->outlets));
      outlet_density_avg = (double *)malloc(sizeof(double) * max(1, lbm->outlets));
      outlet_density_amp = (double *)malloc(sizeof(double) * max(1, lbm->outlets));
      outlet_density_phs = (double *)malloc(sizeof(double) * max(1, lbm->outlets));
      
      for (n = 0; n < lbm->outlets; n++)
	{
	  fscanf (parameters_file, "%le %le %le\n",
		  &outlet_density_avg[n], &outlet_density_amp[n], &outlet_density_phs[n]);
	  
	  outlet_density_avg[n] = lbmConvertPressureToLatticeUnits (outlet_density_avg[n], lbm) / Cs2;
	  outlet_density_amp[n] = lbmConvertPressureGradToLatticeUnits (outlet_density_amp[n], lbm) / Cs2;
	  outlet_density_phs[n] *= DEG_TO_RAD;
	  
	  //if (is_bench)
	  //  {
	  //    outlet_density_avg[ n ] = 1.0;
	  //    outlet_density_amp[ n ] = 0.0;
	  //    outlet_density_phs[ n ] = 0.0;
	  //  }
	}
      lbm_average_inlet_velocity = (double *)malloc(sizeof(double) * lbm->inlets);
      lbm_peak_inlet_velocity    = (double *)malloc(sizeof(double) * lbm->inlets);
      lbm_inlet_normal           = (double *)malloc(sizeof(double) * 3 * lbm->inlets);
      lbm_inlet_count            = (long int *)malloc(sizeof(long int) * lbm->inlets);
      
      if (feof (parameters_file) == 0)
	{
	  is_inlet_normal_available = 1;
	  
	  for (n = 0; n < lbm->inlets; n++)
	    {
	      fscanf (parameters_file, "%le %le %le\n",
		      &lbm_inlet_normal[3*n], &lbm_inlet_normal[3*n+1], &lbm_inlet_normal[3*n+2]);
	    }
	}
      else
	{
	  is_inlet_normal_available = 0;
	}

      fclose (parameters_file);
      
      par_to_send[ 0 ] = 0.1 + (double)lbm->inlets;
      par_to_send[ 1 ] = 0.1 + (double)lbm->outlets;
      par_to_send[ 2 ] = 0.1 + (double)is_inlet_normal_available;
    }
#ifndef NOMPI
  net->err = MPI_Bcast (par_to_send, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
  if (net->id != 0)
    {
      lbm->inlets               = (int)par_to_send[ 0 ];
      lbm->outlets              = (int)par_to_send[ 1 ];
      is_inlet_normal_available = (int)par_to_send[ 2 ];
      
      inlet_density     = (double *)malloc(sizeof(double) * max(1, lbm->inlets));
      inlet_density_avg = (double *)malloc(sizeof(double) * max(1, lbm->inlets));
      inlet_density_amp = (double *)malloc(sizeof(double) * max(1, lbm->inlets));
      inlet_density_phs = (double *)malloc(sizeof(double) * max(1, lbm->inlets));
      
      outlet_density     = (double *)malloc(sizeof(double) * max(1, lbm->outlets));
      outlet_density_avg = (double *)malloc(sizeof(double) * max(1, lbm->outlets));
      outlet_density_amp = (double *)malloc(sizeof(double) * max(1, lbm->outlets));
      outlet_density_phs = (double *)malloc(sizeof(double) * max(1, lbm->outlets));
      
      lbm_average_inlet_velocity = (double *)malloc(sizeof(double) * lbm->inlets);
      lbm_peak_inlet_velocity    = (double *)malloc(sizeof(double) * lbm->inlets);
      lbm_inlet_normal           = (double *)malloc(sizeof(double) * 3 * lbm->inlets);
      lbm_inlet_count            = (long int *)malloc(sizeof(long int) * lbm->inlets);
    }
  else
    {
      for (n = 0; n < lbm->inlets; n++)
	{
	  par_to_send[ 3*n+0 ] = inlet_density_avg[ n ];
	  par_to_send[ 3*n+1 ] = inlet_density_amp[ n ];
	  par_to_send[ 3*n+2 ] = inlet_density_phs[ n ];
	}
      for (n = 0; n < lbm->outlets; n++)
	{
	  par_to_send[ 3*lbm->inlets + 3*n+0 ] = outlet_density_avg[ n ];
	  par_to_send[ 3*lbm->inlets + 3*n+1 ] = outlet_density_amp[ n ];
	  par_to_send[ 3*lbm->inlets + 3*n+2 ] = outlet_density_phs[ n ];
	}
      if (is_inlet_normal_available)
	{
	  for (n = 0; n < lbm->inlets; n++)
	    {
	      par_to_send[ 3*(lbm->inlets+lbm->outlets) + 3*n+0 ] = lbm_inlet_normal[ 3*n+0 ];
	      par_to_send[ 3*(lbm->inlets+lbm->outlets) + 3*n+1 ] = lbm_inlet_normal[ 3*n+1 ];
	      par_to_send[ 3*(lbm->inlets+lbm->outlets) + 3*n+2 ] = lbm_inlet_normal[ 3*n+2 ];
	    }
	}
    }
#ifndef NOMPI
  net->err = MPI_Bcast (par_to_send, 3*(lbm->inlets+lbm->outlets+lbm->inlets), MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
  if (net->id != 0)
    {
      for (n = 0; n < lbm->inlets; n++)
	{
	  inlet_density_avg[ n ] = par_to_send[ 3*n+0 ];
	  inlet_density_amp[ n ] = par_to_send[ 3*n+1 ];
	  inlet_density_phs[ n ] = par_to_send[ 3*n+2 ];
	}
      for (n = 0; n < lbm->outlets; n++)
	{
	  outlet_density_avg[ n ] = par_to_send[ 3*lbm->inlets + 3*n+0 ];
	  outlet_density_amp[ n ] = par_to_send[ 3*lbm->inlets + 3*n+1 ];
	  outlet_density_phs[ n ] = par_to_send[ 3*lbm->inlets + 3*n+2 ];
	}
      if (is_inlet_normal_available)
	{
	  for (n = 0; n < lbm->inlets; n++)
	    {
	      lbm_inlet_normal[ 3*n+0 ] = par_to_send[ 3*(lbm->inlets+lbm->outlets) + 3*n+0 ];
	      lbm_inlet_normal[ 3*n+1 ] = par_to_send[ 3*(lbm->inlets+lbm->outlets) + 3*n+1 ];
	      lbm_inlet_normal[ 3*n+2 ] = par_to_send[ 3*(lbm->inlets+lbm->outlets) + 3*n+2 ];
	    }
	}
    }
  lbmUpdateBoundaryDensities (0, 0, lbm);
  
  lbm->tau = lbmCalculateTau (lbm);
  
  lbm->viscosity = ((2.0 * lbm->tau - 1.0) / 6.0);
  lbm->omega = -1.0 / lbm->tau;
  
  lbm_stress_par = (1.0 - 1.0 / (2.0 * lbm->tau)) / sqrt(2.0);
}


void lbmWriteConfig (int stability, char *output_file_name, LBM *lbm, Net *net)
{
  // this routine writes the flow field on file.
  // the data are collected from the root processor (0 rank).
  // The format comprises:
  // 0- Flag for simulation stability, 0 or 1
  // 1- Voxel size in physical units (units of m)
  // 2- vertex coords of the minimum bounding box with minimum values (x, y and z values)
  // 3- vertex coords of the minimum bounding box with maximum values (x, y and z values)
  // 4- #voxels within the minimum bounding box along the x, y, z axes (3 values)
  // 5- total number of fluid voxels
  // And then a list of the fluid voxels...
  // for each fluid voxel:
  //   a- the (x, y, z) coordinates in lattice units (3 values)
  //   b- the pressure in physical units (mmHg)
  //   c- (x,y,z) components of the velocity field in physical units (3 values, m/s)
  //   d- the von Mises stress in physical units (Pa)
  
  FILE *system_config = NULL;
  XDR	xdr_system_config;
  
  float *local_flow_field, *gathered_flow_field;
  
  double density;
  double pressure;
  double vx, vy, vz;
  double stress;
  double f_eq[15], f_neq[15];
  float pressure_par, velocity_par, stress_par;
  
  int buffer_size;
  int fluid_sites_max;
  int communication_period, communication_iters;
  int period, iters;
  int par;
  int shrinked_sites_x, shrinked_sites_y, shrinked_sites_z;
  int site_i, site_j, site_k;
  int i, j, k;
  int l, m, n;
  int kk;
  
  short int *local_site_data, *gathered_site_data;
  
  unsigned int my_site_id;
  
  // parameters useful to convert pressure, velocity and stress from
  // lattice to physical units
  pressure_par = PULSATILE_PERIOD / (lbm->period * lbm->voxel_size * lbm->voxel_size);
  pressure_par = BLOOD_DENSITY / (mmHg_TO_PASCAL * pressure_par * pressure_par * lbm->voxel_size * lbm->voxel_size);
  velocity_par = 1.0 / (lbm->voxel_size * ((lbm->tau - 0.5) / 3.0) / (BLOOD_VISCOSITY / BLOOD_DENSITY));
  stress_par = ((lbm->tau - 0.5) / 3.0) / (BLOOD_VISCOSITY / BLOOD_DENSITY);
  stress_par = BLOOD_DENSITY / (stress_par * stress_par * lbm->voxel_size * lbm->voxel_size);
  
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
      shrinked_sites_x = 1 + lbm->site_max_x - lbm->site_min_x;
      shrinked_sites_y = 1 + lbm->site_max_y - lbm->site_min_y;
      shrinked_sites_z = 1 + lbm->site_max_z - lbm->site_min_z;
      
      xdr_double (&xdr_system_config, &lbm->voxel_size);
      xdr_int    (&xdr_system_config, &lbm->site_min_x);
      xdr_int    (&xdr_system_config, &lbm->site_min_y);
      xdr_int    (&xdr_system_config, &lbm->site_min_z);
      xdr_int    (&xdr_system_config, &lbm->site_max_x);
      xdr_int    (&xdr_system_config, &lbm->site_max_y);
      xdr_int    (&xdr_system_config, &lbm->site_max_z);
      xdr_int    (&xdr_system_config, &shrinked_sites_x);
      xdr_int    (&xdr_system_config, &shrinked_sites_y);
      xdr_int    (&xdr_system_config, &shrinked_sites_z);
      xdr_int    (&xdr_system_config, &lbm->total_fluid_sites);
    }
  
  fluid_sites_max = 0;
  
  for (n = 0; n < net->procs; n++)
    {
      fluid_sites_max = max(fluid_sites_max, net->fluid_sites[ n ]);
    }
  // "buffer_size" is the size of the flow field buffer to send to the
  // root processor ("local_flow_field") and that to accommodate the
  // received ones from the non-root processors
  // ("gathered_flow_field").  If "buffer_size" is larger the
  // frequency with which data communication to the root processor is
  // performed becomes lower and viceversa
  buffer_size = min(1000000, fluid_sites_max * net->procs);
  
  communication_period = (int)ceil((double)buffer_size / net->procs);
  
  communication_iters = max(1, (int)ceil((double)fluid_sites_max / communication_period));
  
  local_flow_field    = (float *)malloc(sizeof(float) * MACROSCOPIC_PARS * communication_period);
  gathered_flow_field = (float *)malloc(sizeof(float) * MACROSCOPIC_PARS * communication_period * net->procs);
  
  local_site_data    = (short int *)malloc(sizeof(short int) * 3 * communication_period);
  gathered_site_data = (short int *)malloc(sizeof(short int) * 3 * communication_period * net->procs);
  
  for (period = 0; period < communication_period; period++)
    {
      local_site_data[ period*3 ] = -1;
    }
  iters = 0;
  period = 0;
  
  if (!check_conv)
    {
      par = 0;
    }
  else
    {
      par = 1;
    }
  n = -1;
  

  // The following loops scan over every single macrocell (block). If the block is non-empty, it scans the fluid sites within that block
  // If the site is fluid, it calculates the flow field and then is converted to physical units and stored in a buffer to send 
  // to the root processor

  for (i = 0; i < sites_x; i+=block_size)
    {
      for (j = 0; j < sites_y; j+=block_size)
	{
	  for (k = 0; k < sites_z; k+=block_size)
	    {
	      if (net->proc_block[ ++n ].proc_id == NULL)
		{
		  continue;
		}
	      m = -1;
	      
	      for (site_i = i; site_i < i + block_size; site_i++)
		{
		  for (site_j = j; site_j < j + block_size; site_j++)
		    {
		      for (site_k = k; site_k < k + block_size; site_k++)
			{
			  if (net->proc_block[ n ].proc_id[ ++m ] != net->id) continue;
			  
			  my_site_id = net->map_block[ n ].site_data[ m ];
			  
			  if (my_site_id & (1U << 31U)) continue;
			  
			  if (net_site_data[ my_site_id ] == FLUID_TYPE)
			    {
			      lbmFeq (&f_old[ (my_site_id*(par+1)+par)*15 ], &density, &vx, &vy, &vz, f_eq);
			      
			      for (l = 0; l < 15; l++)
				{
				  f_neq[ l ] = f_old[ (my_site_id*(par+1)+par)*15+l ] - f_eq[ l ];
				}
			    }
			  else
			    {
			      lbmCalculateBC (&f_old[ (my_site_id*(par+1)+par)*15 ], net_site_data[ my_site_id ],
					      &density, &vx, &vy, &vz, f_neq);
			    }
			  lbmStress (f_neq, &stress);
			  
			  vx /= density;
			  vy /= density;
			  vz /= density;
			  
			  // conversion from lattice to physical units
			  pressure = REFERENCE_PRESSURE + ((density - 1.0) * Cs2) * pressure_par;
			  vx *= velocity_par;
			  vy *= velocity_par;
			  vz *= velocity_par;
			  stress *= stress_par;
			  
			  local_flow_field[ MACROSCOPIC_PARS*period+0 ] = (float)pressure;
			  local_flow_field[ MACROSCOPIC_PARS*period+1 ] = (float)vx;
			  local_flow_field[ MACROSCOPIC_PARS*period+2 ] = (float)vy;
			  local_flow_field[ MACROSCOPIC_PARS*period+3 ] = (float)vz;
			  local_flow_field[ MACROSCOPIC_PARS*period+4 ] = (float)stress;
			  
			  local_site_data[ period*3+0 ] = site_i;
			  local_site_data[ period*3+1 ] = site_j;
			  local_site_data[ period*3+2 ] = site_k;
			  
			  if (++period != communication_period) continue;
			  
			  period = 0;
			  ++iters;
#ifndef NOMPI
			  net->err = MPI_Gather (local_flow_field, MACROSCOPIC_PARS * communication_period, MPI_FLOAT,
						 gathered_flow_field, MACROSCOPIC_PARS * communication_period, MPI_FLOAT,
						 0, MPI_COMM_WORLD);
			  
			  net->err = MPI_Gather (local_site_data, 3 * communication_period, MPI_SHORT,
						 gathered_site_data, 3 * communication_period, MPI_SHORT, 0, MPI_COMM_WORLD);
#endif
			  if (net->id == 0)
			    {
			      for (l = 0; l < net->procs * communication_period; l++)
				{
				  if (gathered_site_data[ l*3+0 ] == -1) continue;
				  
				  gathered_site_data[ l*3+0 ] -= lbm->site_min_x;
				  gathered_site_data[ l*3+1 ] -= lbm->site_min_y;
				  gathered_site_data[ l*3+2 ] -= lbm->site_min_z;
				  
				  xdr_short (&xdr_system_config, &gathered_site_data[ l*3+0 ]);
				  xdr_short (&xdr_system_config, &gathered_site_data[ l*3+1 ]);
				  xdr_short (&xdr_system_config, &gathered_site_data[ l*3+2 ]);
				  
				  for (kk = 0; kk < MACROSCOPIC_PARS; kk++)
				    {
				      xdr_float (&xdr_system_config, &gathered_flow_field[ MACROSCOPIC_PARS*l+kk ]);
				    }
				}
			    }
			  for (l = 0; l < communication_period; l++)
			    {
			      local_site_data[ l*3 ] = -1;
			    }
			}
		    }
		}
	    }
	}
    }

  if (iters != communication_iters)
    {
      ++iters;
      
      for (; iters <= communication_iters; iters++)
	{
#ifndef NOMPI
	  net->err = MPI_Gather (local_flow_field, MACROSCOPIC_PARS * communication_period, MPI_FLOAT,
				 gathered_flow_field, MACROSCOPIC_PARS * communication_period, MPI_FLOAT,
				 0, MPI_COMM_WORLD);
	  
	  net->err = MPI_Gather (local_site_data, 3 * communication_period, MPI_SHORT,
				 gathered_site_data, 3 * communication_period, MPI_SHORT, 0, MPI_COMM_WORLD);
#endif
      	  if (net->id == 0)
	    {
	      for (l = 0; l < net->procs * communication_period; l++)
		{
		  if (gathered_site_data[ l*3+0 ] == -1) continue;
		  
		  gathered_site_data[ l*3+0 ] -= lbm->site_min_x;
		  gathered_site_data[ l*3+1 ] -= lbm->site_min_y;
		  gathered_site_data[ l*3+2 ] -= lbm->site_min_z;
		  
		  xdr_short (&xdr_system_config, &gathered_site_data[ l*3+0 ]);
		  xdr_short (&xdr_system_config, &gathered_site_data[ l*3+1 ]);
		  xdr_short (&xdr_system_config, &gathered_site_data[ l*3+2 ]);
		  
		  for (kk = 0; kk < MACROSCOPIC_PARS; kk++)
		    {
		      xdr_float (&xdr_system_config, &gathered_flow_field[ MACROSCOPIC_PARS*l+kk ]);
		    }
		}
	    }
	  for (l = 0; l < communication_period; l++)
	    {
	      local_site_data[ l*3 ] = -1;
	    }
	}
    }
  if (net->id == 0)
    {
      fclose (system_config);
    }
  free(gathered_site_data);
  free(local_site_data);
  
  free(gathered_flow_field);
  free(local_flow_field);
}


void lbmWriteConfigASCII (int stability, char *output_file_name, LBM *lbm, Net *net)
{
  // this routine writes the flow field on file.
  // the data are collected from the root processor (0 rank).
  // The format comprises:
  // 0- Flag for simulation stability, 0 or 1
  // 1- Voxel size in physical units (units of m)
  // 2- vertex coords of the minimum bounding box with minimum values (x, y and z values)
  // 3- vertex coords of the minimum bounding box with maximum values (x, y and z values)
  // 4- #voxels within the minimum bounding box along the x, y, z axes (3 values)
  // 5- total number of fluid voxels
  // And then a list of the fluid voxels...
  // for each fluid voxel:
  //   a- the (x, y, z) coordinates in lattice units (3 values)
  //   b- the pressure in physical units (mmHg)
  //   c- (x,y,z) components of the velocity field in physical units (3 values, m/s)
  //   d- the von Mises stress in physical units (Pa)
  
  FILE *system_config = NULL;
  
  float *local_flow_field, *gathered_flow_field;
  
  double density;
  double pressure;
  double vx, vy, vz;
  double stress;
  double f_eq[15], f_neq[15];
  float pressure_par, velocity_par, stress_par;
  
  int buffer_size;
  int fluid_sites_max;
  int communication_period, communication_iters;
  int period, iters;
  int par;
  int shrinked_sites_x, shrinked_sites_y, shrinked_sites_z;
  int site_i, site_j, site_k;
  int i, j, k;
  int l, m, n;
  int kk;
  
  short int *local_site_data, *gathered_site_data;
  
  unsigned int my_site_id;
  
  // parameters useful to convert pressure, velocity and stress from
  // lattice to physical units
  pressure_par = PULSATILE_PERIOD / (lbm->period * lbm->voxel_size * lbm->voxel_size);
  pressure_par = BLOOD_DENSITY / (mmHg_TO_PASCAL * pressure_par * pressure_par * lbm->voxel_size * lbm->voxel_size);
  velocity_par = 1.0 / (lbm->voxel_size * ((lbm->tau - 0.5) / 3.) / (BLOOD_VISCOSITY / BLOOD_DENSITY));
  stress_par = ((lbm->tau - 0.5) / 3.0) / (BLOOD_VISCOSITY / BLOOD_DENSITY);
  stress_par = BLOOD_DENSITY / (stress_par * stress_par * lbm->voxel_size * lbm->voxel_size);
  
  if (net->id == 0)
    {
      system_config = fopen (output_file_name, "w");
      fprintf(system_config, "%i\n", stability);
    }
  
  if (stability == UNSTABLE)
    {
      if (net->id == 0)
	{
	  fclose (system_config);
	}
      return;
    }
  
  if (net->id == 0)
    {
      shrinked_sites_x = 1 + lbm->site_max_x - lbm->site_min_x;
      shrinked_sites_y = 1 + lbm->site_max_y - lbm->site_min_y;
      shrinked_sites_z = 1 + lbm->site_max_z - lbm->site_min_z;
      
      fprintf(system_config, "%e\n", lbm->voxel_size);
      fprintf(system_config, "%i %i %i\n", lbm->site_min_x, lbm->site_min_y, lbm->site_min_z);
      fprintf(system_config, "%i %i %i\n", lbm->site_max_x, lbm->site_max_y, lbm->site_max_z);
      fprintf(system_config, "%i %i %i\n", shrinked_sites_x, shrinked_sites_y, shrinked_sites_z);
      fprintf(system_config, "%i\n", lbm->total_fluid_sites);
    }

  fluid_sites_max = 0;
  
  for (n = 0; n < net->procs; n++)
    {
      fluid_sites_max = max(fluid_sites_max, net->fluid_sites[ n ]);
    }
  // "buffer_size" is the size of the flow field buffer to send to the
  // root processor ("local_flow_field") and that to accommodate the
  // received ones from the non-root processors
  // ("gathered_flow_field").  If "buffer_size" is larger the
  // frequency with which data communication to the root processor is
  // performed becomes lower and viceversa
  buffer_size = min(1000000, fluid_sites_max * net->procs);
  
  communication_period = (int)ceil((double)buffer_size / net->procs);
  
  communication_iters = max(1, (int)ceil((double)fluid_sites_max / communication_period));
  
  local_flow_field    = (float *)malloc(sizeof(float) * MACROSCOPIC_PARS * communication_period);
  gathered_flow_field = (float *)malloc(sizeof(float) * MACROSCOPIC_PARS * communication_period * net->procs);
  
  local_site_data    = (short int *)malloc(sizeof(short int) * 3 * communication_period);
  gathered_site_data = (short int *)malloc(sizeof(short int) * 3 * communication_period * net->procs);
  
  for (period = 0; period < communication_period; period++)
    {
      local_site_data[ period*3 ] = -1;
    }
  iters = 0;
  period = 0;
  
  if (!check_conv)
    {
      par = 0;
    }
  else
    {
      par = 1;
    }
  n = -1;
  

  // The following loops scan over every single macrocell (block). If the block is non-empty, it scans the fluid sites within that block
  // If the site is fluid, it calculates the flow field and then is converted to physical units and stored in a buffer to send 
  // to the root processor

  for (i = 0; i < sites_x; i+=block_size)
    {
      for (j = 0; j < sites_y; j+=block_size)
	{
	  for (k = 0; k < sites_z; k+=block_size)
	    {
	      if (net->proc_block[ ++n ].proc_id == NULL)
		{
		  continue;
		}
	      m = -1;
	      
	      for (site_i = i; site_i < i + block_size; site_i++)
		{
		  for (site_j = j; site_j < j + block_size; site_j++)
		    {
		      for (site_k = k; site_k < k + block_size; site_k++)
			{
			  if (net->proc_block[ n ].proc_id[ ++m ] != net->id) continue;
			  
			  my_site_id = net->map_block[ n ].site_data[ m ];
			  
			  if (my_site_id & (1U << 31U)) continue;
			  
			  if (net_site_data[ my_site_id ] == FLUID_TYPE)
			    {
			      lbmFeq (&f_old[ (my_site_id*(par+1)+par)*15 ], &density, &vx, &vy, &vz, f_eq);
			      
			      for (l = 0; l < 15; l++)
				{
				  f_neq[ l ] = f_old[ (my_site_id*(par+1)+par)*15+l ] - f_eq[ l ];
				}
			    }
			  else
			    {
			      lbmCalculateBC (&f_old[ (my_site_id*(par+1)+par)*15 ], net_site_data[ my_site_id ],
					      &density, &vx, &vy, &vz, f_neq);
			    }
			  lbmStress (f_neq, &stress);
			  
			  vx /= density;
			  vy /= density;
			  vz /= density;
			  
			  // conversion from lattice to physical units
			  pressure = REFERENCE_PRESSURE + ((density - 1.0) * Cs2) * pressure_par;
			  vx *= velocity_par;
			  vy *= velocity_par;
			  vz *= velocity_par;
			  stress *= stress_par;
			  
			  local_flow_field[ MACROSCOPIC_PARS*period+0 ] = (float)pressure;
			  local_flow_field[ MACROSCOPIC_PARS*period+1 ] = (float)vx;
			  local_flow_field[ MACROSCOPIC_PARS*period+2 ] = (float)vy;
			  local_flow_field[ MACROSCOPIC_PARS*period+3 ] = (float)vz;
			  local_flow_field[ MACROSCOPIC_PARS*period+4 ] = (float)stress;
			  
			  local_site_data[ period*3+0 ] = site_i;
			  local_site_data[ period*3+1 ] = site_j;
			  local_site_data[ period*3+2 ] = site_k;
			  
			  if (++period != communication_period) continue;
			  
			  period = 0;
			  ++iters;
#ifndef NOMPI
			  net->err = MPI_Gather (local_flow_field, MACROSCOPIC_PARS * communication_period, MPI_FLOAT,
						 gathered_flow_field, MACROSCOPIC_PARS * communication_period, MPI_FLOAT,
						 0, MPI_COMM_WORLD);
			  
			  net->err = MPI_Gather (local_site_data, 3 * communication_period, MPI_SHORT,
						 gathered_site_data, 3 * communication_period, MPI_SHORT, 0, MPI_COMM_WORLD);
#endif
			  if (net->id == 0)
			    {
			      for (l = 0; l < net->procs * communication_period; l++)
				{
				  if (gathered_site_data[ l*3+0 ] == -1) continue;
				  
				  gathered_site_data[ l*3+0 ] -= lbm->site_min_x;
				  gathered_site_data[ l*3+1 ] -= lbm->site_min_y;
				  gathered_site_data[ l*3+2 ] -= lbm->site_min_z;
				  
				  fprintf(system_config, "%i %i %i",
					  gathered_site_data[ l*3+0 ], gathered_site_data[ l*3+1 ], gathered_site_data[ l*3+2 ]);
				  
				  for (kk = 0; kk < MACROSCOPIC_PARS; kk++)
				    {
			              fprintf(system_config, " %e", gathered_flow_field[ MACROSCOPIC_PARS*l+kk ]);
				    }
				  fprintf(system_config, "\n");
				}
			    }
			  for (l = 0; l < communication_period; l++)
			    {
			      local_site_data[ l*3 ] = -1;
			    }
			}
		    }
		}
	    }
	}
    }

  if (iters != communication_iters)
    {
      ++iters;
      
      for (; iters <= communication_iters; iters++)
	{
#ifndef NOMPI
	  net->err = MPI_Gather (local_flow_field, MACROSCOPIC_PARS * communication_period, MPI_FLOAT,
				 gathered_flow_field, MACROSCOPIC_PARS * communication_period, MPI_FLOAT,
				 0, MPI_COMM_WORLD);
	  
	  net->err = MPI_Gather (local_site_data, 3 * communication_period, MPI_SHORT,
				 gathered_site_data, 3 * communication_period, MPI_SHORT, 0, MPI_COMM_WORLD);
#endif
      	  if (net->id == 0)
	    {
	      for (l = 0; l < net->procs * communication_period; l++)
		{
		  if (gathered_site_data[ l*3+0 ] == -1) continue;
		  
		  gathered_site_data[ l*3+0 ] -= lbm->site_min_x;
		  gathered_site_data[ l*3+1 ] -= lbm->site_min_y;
		  gathered_site_data[ l*3+2 ] -= lbm->site_min_z;
		  
		  fprintf(system_config, "%i %i %i ",
			  gathered_site_data[ l*3+0 ], gathered_site_data[ l*3+1 ], gathered_site_data[ l*3+2 ]);
		  
		  for (kk = 0; kk < MACROSCOPIC_PARS; kk++)
		    {
		      fprintf(system_config, " %e", gathered_flow_field[ MACROSCOPIC_PARS*l+kk ]);
		    }
		  fprintf(system_config, "\n");
		}
	    }
	  for (l = 0; l < communication_period; l++)
	    {
	      local_site_data[ l*3 ] = -1;
	    }
	}
    }
  if (net->id == 0)
    {
      fclose (system_config);
    }
  free(gathered_site_data);
  free(local_site_data);
  
  free(gathered_flow_field);
  free(local_flow_field);
}

