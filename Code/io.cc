/*! \file config.cc
 \brief In this file, the functions useful for the input/output are reported
*/

#include "lb.h"
#include "net.h"
#include "utilityFunctions.h"
#include "xdrReader.h"
#include "xdrFileWriter.h"

//TODO UGH delete when possible
#include "config.h"

#include <limits.h>
#include <sstream>
#include <math.h>
#include <string.h>

using namespace std;

/*!
this function reads the XDR configuration file but does not store the system
and calculate some parameters
*/

void LBM::lbmReadConfig (Net *net) {
  /* Read the config file written by the segtool.
   *
   * All values encoded using XDR format. Uses int, double and u_int.
   * 
   * System parameters:
   *   double stress_type
   *   int blocks_x
   *   int blocks_y
   *   int blocks_z
   *   int block_size 
   *
   * For each block (all blocks_x * blocks_y * blocks_z of them): 
   *
   *   int flag (indicates presence of non-solid sites in the block)
   *   
   *   If flag == 0 go to next block
   *
   *   Otherwise for each site in the block (all block_size^3):
   *
   *     u_int site_data -- this is a bit field which indicates site
   *     type (OR with SITE_TYPE_MASK to get bits zero and one; 00 =
   *     solid, 01 = fluid, 10 = inlet, 11 = outlet) or edgeness (set
   *     bit with PRESSURE_EDGE_MASK)
   * 
   *     If solid or simple fluid, go to next site
   *     
   *     If inlet or outlet (irrespective of edge state) {
   *       double boundary_normal[3]
   *       double boundary_dist
   *     }
   *
   *     If edge bit set {
   *       double wall_normal[3]
   *       double wall_dist
   *     }
   *     
   *     double cut_distances[14]
   */

  FILE* xdrFile = fopen(system_file_name, "r");

  if (xdrFile == NULL) {
    fprintf(stderr, "Unable to open file %s [rank %i], exiting\n", system_file_name, net->id);
    fflush(0x0);
    exit(0x0);
  } else {
    fprintf(stderr, "Opened config file %s [rank %i]\n", system_file_name , net->id);
  }
  fflush(NULL);

  XdrReader myReader = XdrReader(xdrFile);
  
  int i, j, k, ii, jj, kk, l, m, n;
  int flag;
  
  unsigned int site_i, site_j, site_k;
  unsigned int *site_type;
  
  myReader.readDouble(lbm_stress_type);
  myReader.readInt(blocks_x);
  myReader.readInt(blocks_y);
  myReader.readInt(blocks_z);
  myReader.readInt(block_size);
  
  sites_x = blocks_x * block_size;
  sites_y = blocks_y * block_size;
  sites_z = blocks_z * block_size;
  
  sites_in_a_block = block_size * block_size * block_size;

  // shift = log_2(block_size)
  i = block_size;  
  shift = 0;
  while (i > 1) {
    i >>= 1;
    ++shift;
  }
  
  blocks = blocks_x * blocks_y * blocks_z;
  
  net->data_block = (DataBlock *)malloc(sizeof(DataBlock) * blocks);
  
  net->proc_block = (ProcBlock *)malloc(sizeof(ProcBlock) * blocks);
  
  if (lbm_stress_type == SHEAR_STRESS) {
    net->wall_block = (WallBlock *)malloc(sizeof(WallBlock) * blocks);
  }
  total_fluid_sites = 0;
  
  site_min_x = INT_MAX;
  site_min_y = INT_MAX;
  site_min_z = INT_MAX;
  site_max_x = INT_MIN;
  site_max_y = INT_MIN;
  site_max_z = INT_MIN;
  
  net->fr_time = UtilityFunctions::myClock ();
  
  n = -1;
  
  for (i = 0; i < blocks_x; i++) {
    for (j = 0; j < blocks_y; j++) {
      for (k = 0; k < blocks_z; k++) {
	++n;
	
	net->data_block[n].site_data = NULL;
	net->proc_block[n].proc_id   = NULL;
	net->wall_block[n].wall_data = NULL;
	
	myReader.readInt(flag);
	
	if (flag == 0) continue;
	// Block contains some non-solid sites
	
	net->data_block[n].site_data = (unsigned int *)malloc(sizeof(unsigned int) * sites_in_a_block);
	net->proc_block[n].proc_id   = (int *)malloc(sizeof(int) * sites_in_a_block);
	
	m = -1;
	
	for (ii = 0; ii < block_size; ii++) {
	  site_i = (i << shift) + ii;
	  
	  for (jj = 0; jj < block_size; jj++) {
	    site_j = (j << shift) + jj;
	    
	    for (kk = 0; kk < block_size; kk++) {
	      site_k = (k << shift) + kk;
	      
	      ++m;
	      
	      site_type =  &net->data_block[n].site_data[m];
	      myReader.readUnsignedInt(*site_type);
	      
	      if ((*site_type & SITE_TYPE_MASK) == SOLID_TYPE) {
		net->proc_block[n].proc_id[m] = 1 << 30;
		continue;
	      }
	      net->proc_block[n].proc_id[ m ] = -1;
	      
	      ++total_fluid_sites;
	      
	      site_min_x = UtilityFunctions::min(site_min_x, site_i);
	      site_min_y = UtilityFunctions::min(site_min_y, site_j);
	      site_min_z = UtilityFunctions::min(site_min_z, site_k);
	      site_max_x = UtilityFunctions::max(site_max_x, site_i);
	      site_max_y = UtilityFunctions::max(site_max_y, site_j);
	      site_max_z = UtilityFunctions::max(site_max_z, site_k);
	      
	      if (lbm_stress_type == SHEAR_STRESS &&
		  lbmCollisionType (*site_type) != FLUID) {
		// Neither solid nor simple fluid
		if (net->wall_block[n].wall_data == NULL) {
		  net->wall_block[n].wall_data = (WallData *)malloc(sizeof(WallData) * sites_in_a_block);
		}
		
		if (lbmCollisionType (*site_type) & INLET ||
		    lbmCollisionType (*site_type) & OUTLET) {
		  // INLET or OUTLET or both
		  for (l = 0; l < 3; l++)
		    myReader.readDouble(net->wall_block[n].wall_data[m].boundary_nor[l]);
		  
		  myReader.readDouble(net->wall_block[n].wall_data[m].boundary_dist);
		}
		
		if (lbmCollisionType(*site_type) & EDGE) {
		  // EDGE bit set
		  for (l = 0; l < 3; l++)
		    myReader.readDouble(net->wall_block[n].wall_data[m].wall_nor[l]);
		  
		  myReader.readDouble(net->wall_block[n].wall_data[m].wall_dist);
		}
		
		for (l = 0; l < 14; l++)
		  myReader.readDouble(net->wall_block[n].wall_data[m].cut_dist[l]);
	      }
	    } // kk
	  }   // jj
	}     // ii
      } // k
    }   // j
  }     // i
  
  fclose (xdrFile);

  net->fr_time = UtilityFunctions::myClock () - net->fr_time;
}

/*!
through this function the processor 0 reads the LB parameters
and then communicate them to the other processors
*/
void LBM::lbmReadParameters (char *parameters_file_name, Net *net)
{
  
  double par_to_send[10000];
  double nor[3], pos[3];
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
      
      fscanf (parameters_file, "%i\n", &inlets);

      lbmInitialiseInlets(inlets);

      for (n = 0; n < inlets; n++)
	{
	  fscanf (parameters_file, "%le %le %le\n",
		  &inlet_density_avg[n], &inlet_density_amp[n], &inlet_density_phs[n]);
	  
	  inlet_density_avg[n] = lbmConvertPressureToLatticeUnits (inlet_density_avg[n]) / Cs2;
	  inlet_density_amp[n] = lbmConvertPressureGradToLatticeUnits (inlet_density_amp[n]) / Cs2;
	  inlet_density_phs[n] *= DEG_TO_RAD;
	}
      fscanf (parameters_file, "%i\n", &outlets);


      lbmInitialiseOutlets(outlets);
      
      for (n = 0; n < outlets; n++)
	{
	  fscanf (parameters_file, "%le %le %le\n",
		  &outlet_density_avg[n], &outlet_density_amp[n], &outlet_density_phs[n]);
	  
	  outlet_density_avg[n] = lbmConvertPressureToLatticeUnits (outlet_density_avg[n]) / Cs2;
	  outlet_density_amp[n] = lbmConvertPressureGradToLatticeUnits (outlet_density_amp[n]) / Cs2;
	  outlet_density_phs[n] *= DEG_TO_RAD;
	}
      lbm_average_inlet_velocity = (double *)malloc(sizeof(double) * inlets);
      lbm_peak_inlet_velocity    = (double *)malloc(sizeof(double) * inlets);
      lbm_inlet_normal           = (double *)malloc(sizeof(double) * 3 * inlets);
      lbm_inlet_count            = (long int *)malloc(sizeof(long int) * inlets);
      
      if (feof (parameters_file) == 0)
	{
	  is_inlet_normal_available = 1;
	  
	  for (n = 0; n < inlets; n++)
	    fscanf (parameters_file, "%le %le %le\n",
		      &lbm_inlet_normal[3*n], &lbm_inlet_normal[3*n+1], &lbm_inlet_normal[3*n+2]);
	}
      else
	{
	  is_inlet_normal_available = 0;
	}
      if (feof (parameters_file) == 0)
	{
	  for (n = 0; n < outlets; n++)
	    fscanf (parameters_file, "%le %le %le\n", &nor[0], &nor[1], &nor[2]);
	}
      if (feof (parameters_file) == 0)
	{
	  for (n = 0; n < inlets; n++)
	    fscanf (parameters_file, "%le %le %le\n", &pos[0], &pos[1], &pos[2]);
	}
      if (feof (parameters_file) == 0)
	{
	  for (n = 0; n < outlets; n++)
	    fscanf (parameters_file, "%le %le %le\n", &pos[0], &pos[1], &pos[2]);
	}
      fclose (parameters_file);
      
      par_to_send[ 0 ] = 0.1 + (double)inlets;
      par_to_send[ 1 ] = 0.1 + (double)outlets;
      par_to_send[ 2 ] = 0.1 + (double)is_inlet_normal_available;
    }
#ifndef NOMPI
  net->err = MPI_Bcast (par_to_send, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
  if (net->id != 0)
    {
      inlets               = (int)par_to_send[ 0 ];
      outlets              = (int)par_to_send[ 1 ];
      is_inlet_normal_available = (int)par_to_send[ 2 ];

      lbmInitialiseInlets(inlets); 
      lbmInitialiseOutlets(outlets);
      
      lbm_average_inlet_velocity = (double *)malloc(sizeof(double) * inlets);
      lbm_peak_inlet_velocity    = (double *)malloc(sizeof(double) * inlets);
      lbm_inlet_normal           = (double *)malloc(sizeof(double) * 3 * inlets);
      lbm_inlet_count            = (long int *)malloc(sizeof(long int) * inlets);
    }
  else
    {
      for (n = 0; n < inlets; n++)
	{
	  par_to_send[ 3*n+0 ] = inlet_density_avg[ n ];
	  par_to_send[ 3*n+1 ] = inlet_density_amp[ n ];
	  par_to_send[ 3*n+2 ] = inlet_density_phs[ n ];
	}
      for (n = 0; n < outlets; n++)
	{
	  par_to_send[ 3*inlets + 3*n+0 ] = outlet_density_avg[ n ];
	  par_to_send[ 3*inlets + 3*n+1 ] = outlet_density_amp[ n ];
	  par_to_send[ 3*inlets + 3*n+2 ] = outlet_density_phs[ n ];
	}
      if (is_inlet_normal_available)
	{
	  for (n = 0; n < inlets; n++)
	    {
	      par_to_send[ 3*(inlets+outlets) + 3*n+0 ] = lbm_inlet_normal[ 3*n+0 ];
	      par_to_send[ 3*(inlets+outlets) + 3*n+1 ] = lbm_inlet_normal[ 3*n+1 ];
	      par_to_send[ 3*(inlets+outlets) + 3*n+2 ] = lbm_inlet_normal[ 3*n+2 ];
	    }
	}
    }
#ifndef NOMPI
  net->err = MPI_Bcast (par_to_send, 3*(inlets+outlets+inlets), MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
  if (net->id != 0)
    {
      for (n = 0; n < inlets; n++)
	{
	  inlet_density_avg[ n ] = par_to_send[ 3*n+0 ];
	  inlet_density_amp[ n ] = par_to_send[ 3*n+1 ];
	  inlet_density_phs[ n ] = par_to_send[ 3*n+2 ];
	}
      for (n = 0; n < outlets; n++)
	{
	  outlet_density_avg[ n ] = par_to_send[ 3*inlets + 3*n+0 ];
	  outlet_density_amp[ n ] = par_to_send[ 3*inlets + 3*n+1 ];
	  outlet_density_phs[ n ] = par_to_send[ 3*inlets + 3*n+2 ];
	}
      if (is_inlet_normal_available)
	{
	  for (n = 0; n < inlets; n++)
	    {
	      lbm_inlet_normal[ 3*n+0 ] = par_to_send[ 3*(inlets+outlets) + 3*n+0 ];
	      lbm_inlet_normal[ 3*n+1 ] = par_to_send[ 3*(inlets+outlets) + 3*n+1 ];
	      lbm_inlet_normal[ 3*n+2 ] = par_to_send[ 3*(inlets+outlets) + 3*n+2 ];
	    }
	}
    }
  lbmUpdateBoundaryDensities (0, 0);
  
  RecalculateTauViscosityOmega ();
}

void LBM::lbmInitialiseInlets(int numberOfInlets)
{     
  inlet_density     = (double *)malloc(sizeof(double) * UtilityFunctions::max(1, numberOfInlets));
  inlet_density_avg = (double *)malloc(sizeof(double) * UtilityFunctions::max(1, numberOfInlets));
  inlet_density_amp = (double *)malloc(sizeof(double) * UtilityFunctions::max(1, numberOfInlets));
  inlet_density_phs = (double *)malloc(sizeof(double) * UtilityFunctions::max(1, numberOfInlets));
}

void LBM::lbmInitialiseOutlets(int numberOfOutlets)
{
  outlet_density     = (double *)malloc(sizeof(double) * UtilityFunctions::max(1, numberOfOutlets));
  outlet_density_avg = (double *)malloc(sizeof(double) * UtilityFunctions::max(1, numberOfOutlets));
  outlet_density_amp = (double *)malloc(sizeof(double) * UtilityFunctions::max(1, numberOfOutlets));
  outlet_density_phs = (double *)malloc(sizeof(double) * UtilityFunctions::max(1, numberOfOutlets));
}

void LBM::lbmWriteConfig (int stability, char *output_file_name, Net *net)
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
  //   d- the von Mises stress or shear stress in physical units (Pa)
  //      (the stored shear stress is equal to -1 if the fluid voxel is not at the wall)
  
  XdrWriter* myWriter;
  
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
  pressure_par = PULSATILE_PERIOD / (period * voxel_size * voxel_size);
  pressure_par = BLOOD_DENSITY / (mmHg_TO_PASCAL * pressure_par * pressure_par * voxel_size * voxel_size);
  velocity_par = 1.0 / (voxel_size * ((tau - 0.5) / 3.0) / (BLOOD_VISCOSITY / BLOOD_DENSITY));
  stress_par = ((tau - 0.5) / 3.0) / (BLOOD_VISCOSITY / BLOOD_DENSITY);
  stress_par = BLOOD_DENSITY / (stress_par * stress_par * voxel_size * voxel_size);
  
  if (net->id == 0)
    {
      myWriter = new XdrFileWriter(output_file_name);

      myWriter->writeInt(&stability);
  
      if (stability == UNSTABLE)
      {
        return;
      }
  
      shrinked_sites_x = 1 + site_max_x - site_min_x;
      shrinked_sites_y = 1 + site_max_y - site_min_y;
      shrinked_sites_z = 1 + site_max_z - site_min_z;
      
      myWriter->writeDouble(&voxel_size);
      myWriter->writeInt(&site_min_x);
      myWriter->writeInt(&site_min_y);
      myWriter->writeInt(&site_min_z);
      myWriter->writeInt(&site_max_x);
      myWriter->writeInt(&site_max_y);
      myWriter->writeInt(&site_max_z);
      myWriter->writeInt(&shrinked_sites_x);
      myWriter->writeInt(&shrinked_sites_y);
      myWriter->writeInt(&shrinked_sites_z);
      myWriter->writeInt(&total_fluid_sites);
    }
  
  fluid_sites_max = 0;
  
  for (n = 0; n < net->procs; n++)
    {
      fluid_sites_max = UtilityFunctions::max(fluid_sites_max, net->fluid_sites[ n ]);
    }
  // "buffer_size" is the size of the flow field buffer to send to the
  // root processor ("local_flow_field") and that to accommodate the
  // received ones from the non-root processors
  // ("gathered_flow_field").  If "buffer_size" is larger the
  // frequency with which data communication to the root processor is
  // performed becomes lower and viceversa
  buffer_size = UtilityFunctions::min(1000000, fluid_sites_max * net->procs);
  
  communication_period = (int)ceil((double)buffer_size / net->procs);
  
  communication_iters = UtilityFunctions::max(1, (int)ceil((double)fluid_sites_max / communication_period));
  
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
			  
			  if (net->net_site_data[ my_site_id ] == FLUID_TYPE)
			    {
			      lbmFeq (&f_old[ (my_site_id*(par+1)+par)*15 ], &density, &vx, &vy, &vz, f_eq);
			      
			      for (l = 0; l < 15; l++)
				{
				  f_neq[ l ] = f_old[ (my_site_id*(par+1)+par)*15+l ] - f_eq[ l ];
				}
			    }
			  else
			    {
			      lbmCalculateBC (&f_old[ (my_site_id*(par+1)+par)*15 ], net->net_site_data[ my_site_id ],
					      &density, &vx, &vy, &vz, f_neq);
			    }
			  if (lbm_stress_type == SHEAR_STRESS)
			    {
			      if (net->net_site_nor[ my_site_id*3 ] >= 1.0e+30)
				{
				  stress = -1.0;
				}
			      else
				{
				  lbmStress (density, f_neq, &net->net_site_nor[ my_site_id*3 ], &stress);
				}
			    }
			  else
			    {
			      lbmStress (f_neq, &stress);
			    }
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
				  
				  gathered_site_data[ l*3+0 ] -= site_min_x;
				  gathered_site_data[ l*3+1 ] -= site_min_y;
				  gathered_site_data[ l*3+2 ] -= site_min_z;
				  
				  myWriter->writeShort(&gathered_site_data[ l*3+0 ]);
				  myWriter->writeShort(&gathered_site_data[ l*3+1 ]);
				  myWriter->writeShort(&gathered_site_data[ l*3+2 ]);
				  
				  for (kk = 0; kk < MACROSCOPIC_PARS; kk++)
				    {
				      myWriter->writeFloat(&gathered_flow_field[ MACROSCOPIC_PARS*l+kk ]);
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
		  
		  gathered_site_data[ l*3+0 ] -= site_min_x;
		  gathered_site_data[ l*3+1 ] -= site_min_y;
		  gathered_site_data[ l*3+2 ] -= site_min_z;
		  
		  myWriter->writeShort(&gathered_site_data[ l*3+0 ]);
		  myWriter->writeShort(&gathered_site_data[ l*3+1 ]);
		  myWriter->writeShort(&gathered_site_data[ l*3+2 ]);
		  
		  for (kk = 0; kk < MACROSCOPIC_PARS; kk++)
		    {
		      myWriter->writeFloat(&gathered_flow_field[ MACROSCOPIC_PARS*l+kk ]);
		    }
		}
	    }
	  for (l = 0; l < communication_period; l++)
	    {
	      local_site_data[ l*3 ] = -1;
	    }
	}
    }

  free(gathered_site_data);
  free(local_site_data);
  
  free(gathered_flow_field);
  free(local_flow_field);

  delete myWriter;
}


void LBM::lbmWriteConfigASCII (int stability, char *output_file_name, Net *net)
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
  //      (the stored shear stress is equal to -1 if the fluid voxel is not at the wall)
  
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
  int comPeriodDelta, iters;
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
  pressure_par = PULSATILE_PERIOD / (period * voxel_size * voxel_size);
  pressure_par = BLOOD_DENSITY / (mmHg_TO_PASCAL * pressure_par * pressure_par * voxel_size * voxel_size);
  velocity_par = 1.0 / (voxel_size * ((tau - 0.5) / 3.) / (BLOOD_VISCOSITY / BLOOD_DENSITY));
  stress_par = ((tau - 0.5) / 3.0) / (BLOOD_VISCOSITY / BLOOD_DENSITY);
  stress_par = BLOOD_DENSITY / (stress_par * stress_par * voxel_size * voxel_size);
  
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
      shrinked_sites_x = 1 + site_max_x - site_min_x;
      shrinked_sites_y = 1 + site_max_y - site_min_y;
      shrinked_sites_z = 1 + site_max_z - site_min_z;
      
      fprintf(system_config, "%e\n", voxel_size);
      fprintf(system_config, "%i %i %i\n", site_min_x, site_min_y, site_min_z);
      fprintf(system_config, "%i %i %i\n", site_max_x, site_max_y, site_max_z);
      fprintf(system_config, "%i %i %i\n", shrinked_sites_x, shrinked_sites_y, shrinked_sites_z);
      fprintf(system_config, "%i\n", total_fluid_sites);
    }

  fluid_sites_max = 0;
  
  for (n = 0; n < net->procs; n++)
    {
      fluid_sites_max = UtilityFunctions::max(fluid_sites_max, net->fluid_sites[ n ]);
    }
  // "buffer_size" is the size of the flow field buffer to send to the
  // root processor ("local_flow_field") and that to accommodate the
  // received ones from the non-root processors
  // ("gathered_flow_field").  If "buffer_size" is larger the
  // frequency with which data communication to the root processor is
  // performed becomes lower and viceversa
  buffer_size = UtilityFunctions::min(1000000, fluid_sites_max * net->procs);
  
  communication_period = (int)ceil((double)buffer_size / net->procs);
  
  communication_iters = UtilityFunctions::max(1, (int)ceil((double)fluid_sites_max / communication_period));
  
  local_flow_field    = (float *)malloc(sizeof(float) * MACROSCOPIC_PARS * communication_period);
  gathered_flow_field = (float *)malloc(sizeof(float) * MACROSCOPIC_PARS * communication_period * net->procs);
  
  local_site_data    = (short int *)malloc(sizeof(short int) * 3 * communication_period);
  gathered_site_data = (short int *)malloc(sizeof(short int) * 3 * communication_period * net->procs);
  
  for (comPeriodDelta = 0; comPeriodDelta < communication_period; comPeriodDelta++)
    {
      local_site_data[ comPeriodDelta*3 ] = -1;
    }
  iters = 0;
  comPeriodDelta = 0;
  
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
			  
			  if (net->net_site_data[ my_site_id ] == FLUID_TYPE)
			    {
			      lbmFeq (&f_old[ (my_site_id*(par+1)+par)*15 ], &density, &vx, &vy, &vz, f_eq);
			      
			      for (l = 0; l < 15; l++)
				{
				  f_neq[ l ] = f_old[ (my_site_id*(par+1)+par)*15+l ] - f_eq[ l ];
				}
			    }
			  else
			    {
			      lbmCalculateBC (&f_old[ (my_site_id*(par+1)+par)*15 ], net->net_site_data[ my_site_id ],
					      &density, &vx, &vy, &vz, f_neq);
			    }
			  if (lbm_stress_type == SHEAR_STRESS)
			    {
			      if (net->net_site_nor[ my_site_id*3 ] >= 1.0e+30)
				{
				  stress = -1.0;
				}
			      else
				{
				  lbmStress (density, f_neq, &net->net_site_nor[ my_site_id*3 ], &stress);
				}
			    }
			  else
			    {
			      lbmStress (f_neq, &stress);
			    }
			  vx /= density;
			  vy /= density;
			  vz /= density;
			  
			  // conversion from lattice to physical units
			  pressure = REFERENCE_PRESSURE + ((density - 1.0) * Cs2) * pressure_par;
			  vx *= velocity_par;
			  vy *= velocity_par;
			  vz *= velocity_par;
			  stress *= stress_par;
			  
			  local_flow_field[ MACROSCOPIC_PARS*comPeriodDelta+0 ] = (float)pressure;
			  local_flow_field[ MACROSCOPIC_PARS*comPeriodDelta+1 ] = (float)vx;
			  local_flow_field[ MACROSCOPIC_PARS*comPeriodDelta+2 ] = (float)vy;
			  local_flow_field[ MACROSCOPIC_PARS*comPeriodDelta+3 ] = (float)vz;
			  local_flow_field[ MACROSCOPIC_PARS*comPeriodDelta+4 ] = (float)stress;
			  
			  local_site_data[ comPeriodDelta*3+0 ] = site_i;
			  local_site_data[ comPeriodDelta*3+1 ] = site_j;
			  local_site_data[ comPeriodDelta*3+2 ] = site_k;
			  
			  if (++comPeriodDelta != communication_period) continue;
			  
			  comPeriodDelta = 0;
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
				  
				  gathered_site_data[ l*3+0 ] -= site_min_x;
				  gathered_site_data[ l*3+1 ] -= site_min_y;
				  gathered_site_data[ l*3+2 ] -= site_min_z;
				  
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
		  
		  gathered_site_data[ l*3+0 ] -= site_min_x;
		  gathered_site_data[ l*3+1 ] -= site_min_y;
		  gathered_site_data[ l*3+2 ] -= site_min_z;
		  
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

