#include "streaklineDrawer.h"
// TODO this could probably be reduced to the net class and some visualisation class.
#include "rt.h"
#include "utilityFunctions.h"
// TODO: Mixture of malloc and new is a bad smell, and just using new / delete would be much nicer.
#include <stdlib.h>
#include <math.h>

void streaklineDrawer::slInitializeVelFieldBlock (int site_i, int site_j, int site_k, int proc_id)
{
  if (site_i < 0 || site_i >= sites_x ||
      site_j < 0 || site_j >= sites_y ||
      site_k < 0 || site_k >= sites_z)
    {
      return;
    }
  int i = site_i >> shift;
  int j = site_j >> shift;
  int k = site_k >> shift;
  
  int block_id =  (i * blocks_y + j) * blocks_z + k;
  int site_id;
  
  if (velocity_field[ block_id ].vel_site_data == NULL)
    {
      velocity_field[ block_id ].vel_site_data = new VelSiteData[sites_in_a_block];
      
      for (site_id = 0; site_id < sites_in_a_block; site_id++)
	{
	  velocity_field[ block_id ].vel_site_data[ site_id ].proc_id = -1;
	  velocity_field[ block_id ].vel_site_data[ site_id ].counter = 0;
	}
    }
  int ii = site_i - (i << shift);
  int jj = site_j - (j << shift);
  int kk = site_k - (k << shift);
  
  site_id = (((ii << shift) + jj) << shift) + kk;
  velocity_field[ block_id ].vel_site_data[ site_id ].proc_id = proc_id;
}


streaklineDrawer::VelSiteData *streaklineDrawer::slVelSiteDataPointer (int site_i, int site_j, int site_k)
{
  if (site_i < 0 || site_i >= sites_x ||
      site_j < 0 || site_j >= sites_y ||
      site_k < 0 || site_k >= sites_z)
    {
      return NULL;
    }
  int i = site_i >> shift;
  int j = site_j >> shift;
  int k = site_k >> shift;
  
  int block_id = (i * blocks_y + j) * blocks_z + k;
  
  if (velocity_field[ block_id ].vel_site_data == NULL)
    {
      return NULL;
    }
  int ii = site_i - (i << shift);
  int jj = site_j - (j << shift);
  int kk = site_k - (k << shift);
  
  int site_id = (((ii << shift) + jj) << shift) + kk;
  
  return &velocity_field[ block_id ].vel_site_data[ site_id ];
}


void streaklineDrawer::slParticleVelocity (Particle *particle_p, float v[2][2][2][3], float interp_v[3])
{
  float dx, dy, dz;
  float v_00z, v_01z, v_10z, v_11z, v_0y, v_1y;
  
  
  dx = particle_p->x - (int)particle_p->x;
  dy = particle_p->y - (int)particle_p->y;
  dz = particle_p->z - (int)particle_p->z;
  
  for (int l = 0; l < 3; l++)
    {
      v_00z = (1.F - dz) * v[0][0][0][l] + dz * v[0][0][1][l];
      v_01z = (1.F - dz) * v[0][1][0][l] + dz * v[0][1][1][l];
      v_10z = (1.F - dz) * v[1][0][0][l] + dz * v[1][0][1][l];
      v_11z = (1.F - dz) * v[1][1][0][l] + dz * v[1][1][1][l];
      
      v_0y = (1.F - dy) * v_00z + dy * v_01z;
      v_1y = (1.F - dy) * v_10z + dy * v_11z;
      
      interp_v[l] = (1.F - dx) * v_0y + dx * v_1y;
    }
}


void streaklineDrawer::slCreateParticle (float x, float y, float z, float vel, int inlet_id)
{
  if (particles == particles_max)
    {
      particles_max *= 2;
      particle = (Particle *)realloc(particle, sizeof(Particle) * particles_max);
    }
  particle[ particles ].x        = x;
  particle[ particles ].y        = y;
  particle[ particles ].z        = z;
  particle[ particles ].vel      = vel;
  particle[ particles ].inlet_id = inlet_id;
  ++particles;
}


void streaklineDrawer::slDeleteParticle (int p_index)
{
  if (--particles <= 0) return;
  
  // its data are reslaced with those of the last particle;
  if (p_index != particles)
    {
      particle[ p_index ].x        = particle[ particles ].x;
      particle[ p_index ].y        = particle[ particles ].y;
      particle[ p_index ].z        = particle[ particles ].z;
      particle[ p_index ].vx       = particle[ particles ].vx;
      particle[ p_index ].vy       = particle[ particles ].vy;
      particle[ p_index ].vz       = particle[ particles ].vz;
      particle[ p_index ].vel      = particle[ particles ].vel;
      particle[ p_index ].inlet_id = particle[ particles ].inlet_id;
    }
}


void streaklineDrawer::slCreateSeedParticles ()
{
  for (int n = 0; n < particle_seeds; n++)
    {
      slCreateParticle (particle_seed[n].x,
			particle_seed[n].y,
			particle_seed[n].z,
			0.0F,
			particle_seed[n].inlet_id);
    }
}


void streaklineDrawer::slLocalVelField (int p_index, float v[2][2][2][3], int *is_interior, Net *net)
{
  double density;
  double vx, vy, vz;
  
  int site_i, site_j, site_k;
  int neigh_i, neigh_j, neigh_k;
  int i, j, k;
  int m;
  int c1, c2;
  
  VelSiteData *vel_site_data_p;
  
  
  site_i = (int)particle[ p_index ].x;
  site_j = (int)particle[ p_index ].y;
  site_k = (int)particle[ p_index ].z;
  
  *is_interior = 1;
  
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
  for (i = 0; i < 2; i++)
    {
      neigh_i = site_i + i;
      
      for (j = 0; j < 2; j++)
	{
	  neigh_j = site_j + j;
	  
	  for (k = 0; k < 2; k++)
	    {
	      neigh_k = site_k + k;
	      
	      vel_site_data_p = slVelSiteDataPointer (neigh_i, neigh_j, neigh_k);
	      
	      if (vel_site_data_p == NULL ||
		  vel_site_data_p->proc_id == -1)
		{
		  // it is a solid site and the velocity is
		  // assumed to be zero
		  v[i][j][k][0] = v[i][j][k][1] = v[i][j][k][2] = 0.0F;
		  continue;
		}
	      if (vel_site_data_p->proc_id != net->id)
		{
		  *is_interior = 0;
		}
	      if (vel_site_data_p->counter == counter)
		{
		  // it means that the local velocity has already been
		  // calculated at the current time step if the site
		  // belongs to the current processor; if not, the
		  // following instructions have no effect
		  v[i][j][k][0] = vel_site_data_p->vx;
		  v[i][j][k][1] = vel_site_data_p->vy;
		  v[i][j][k][2] = vel_site_data_p->vz;
		}
	      else if (vel_site_data_p->proc_id == net->id)
		{
		  // the local counter is set equal to the global one
		  // and the local velocity is calculated
		  vel_site_data_p->counter = counter;
		  
		  lbmDensityAndVelocity (&f_old[ vel_site_data_p->site_id*c1+c2 ],
					 &density, &vx, &vy, &vz);
		  
		  v[i][j][k][0] = vel_site_data_p->vx = vx / density;
		  v[i][j][k][1] = vel_site_data_p->vy = vy / density;
		  v[i][j][k][2] = vel_site_data_p->vz = vz / density;
		}
	      else
		{
		  vel_site_data_p->counter = counter;
		  
		  m = from_proc_id_to_neigh_proc_index[ vel_site_data_p->proc_id ];
		  
		  neigh_proc[m].s_to_send[ 3*neigh_proc[m].send_vs+0 ] = neigh_i;
		  neigh_proc[m].s_to_send[ 3*neigh_proc[m].send_vs+1 ] = neigh_j;
		  neigh_proc[m].s_to_send[ 3*neigh_proc[m].send_vs+2 ] = neigh_k;
		  ++neigh_proc[m].send_vs;
		}
	    }
	}
    }
}


void streaklineDrawer::slInit (Net *net)
{
  int site_i, site_j, site_k;
  int neigh_i, neigh_j, neigh_k;
  int i, j, k;
  int m, mm, n;
  int flag;
  int inlet_sites;
  int *neigh_proc_id;
  
  unsigned int site_data;
  
  DataBlock *map_block_p;
  
  ProcBlock *proc_block_p;
  
  VelSiteData *vel_site_data_p;
  
  
  particles_max = 10000;
  particle = new Particle[particles_max];
  particles = 0;
  
  particle_seeds_max = 100;
  particle_seed = new Particle[particle_seeds_max];
  particle_seeds = 0;
  
  neigh_procs = 0;
  
  velocity_field = new VelocityField [blocks];
  
  for (n = 0; n < blocks; n++)
    {
      velocity_field[ n ].vel_site_data = NULL;
    }
  for (m = 0; m <  NEIGHBOUR_PROCS_MAX; m++)
    {
      neigh_proc[ m ].send_vs = 0;
    }
  counter = 1;
  inlet_sites = 0;
  n = -1;
  
  for (i = 0; i < sites_x; i += block_size)
    for (j = 0; j < sites_y; j += block_size)
      for (k = 0; k < sites_z; k += block_size)
	{
	  map_block_p = &net->map_block[ ++n ];
	  
	  if (map_block_p->site_data == NULL) continue;
	  
	  proc_block_p = &net->proc_block[ n ];
	  
	  m = -1;
	  
	  for (site_i = i; site_i < i + block_size; site_i++)
	    for (site_j = j; site_j < j + block_size; site_j++)
	      for (site_k = k; site_k < k + block_size; site_k++)
		{
		  if (proc_block_p->proc_id[ ++m ] != net->id) continue;
		  
		  for (neigh_i = UtilityFunctions::max(0, site_i-1); neigh_i <= UtilityFunctions::min(sites_x-1, site_i+1); neigh_i++)
		    for (neigh_j = UtilityFunctions::max(0, site_j-1); neigh_j <= UtilityFunctions::min(sites_y-1, site_j+1); neigh_j++)
		      for (neigh_k = UtilityFunctions::max(0, site_k-1); neigh_k <= UtilityFunctions::min(sites_z-1, site_k+1); neigh_k++)
			{
			  neigh_proc_id = net->netProcIdPointer (neigh_i, neigh_j, neigh_k);
			  
			  if (neigh_proc_id == NULL || *neigh_proc_id == (1 << 30))
			    {
			      continue;
			    }
			  slInitializeVelFieldBlock (neigh_i, neigh_j, neigh_k, *neigh_proc_id);
			  
			  if (*neigh_proc_id == net->id) continue;
			  
			  vel_site_data_p = slVelSiteDataPointer (neigh_i, neigh_j, neigh_k);
			  
			  if (vel_site_data_p->counter == counter) continue;
			  
			  vel_site_data_p->counter = counter;
			  
			  for (mm = 0, flag = 1; mm < neigh_procs && flag; mm++)
			    {
			      if (*neigh_proc_id == neigh_proc[ mm ].id)
				{
				  flag = 0;
				  ++neigh_proc[ mm ].send_vs;
				}
			    }
			  if (!flag) continue;
			  
			  if (neigh_procs == NEIGHBOUR_PROCS_MAX)
			    {
			      printf (" too many inter processor neighbours in slInit()\n");
			      printf (" the execution is terminated\n");
#ifndef NOMPI
			      net->err = MPI_Abort (MPI_COMM_WORLD, 1);
#else
			      exit(1);
#endif
			    }
			  neigh_proc[ neigh_procs ].id = *neigh_proc_id;
			  neigh_proc[ neigh_procs ].send_vs = 1;
			  ++neigh_procs;
			}
		  site_data = net->net_site_data[ map_block_p->site_data[m] ];
		  
		  // if the lattice site is an not inlet one
		  if ((site_data & SITE_TYPE_MASK) != INLET_TYPE)
		    {
		      continue;
		    }
		  ++inlet_sites;
		  
		  if (inlet_sites%50 != 0) continue;
		  
		  if (particle_seeds == particle_seeds_max)
		    {
		      particle_seeds_max *= 2;
		      particle_seed = (Particle *)realloc(particle_seed,
							      sizeof(Particle) * particle_seeds_max);
		    }
		  particle_seed[ particle_seeds ].x = (float)site_i;
		  particle_seed[ particle_seeds ].y = (float)site_j;
		  particle_seed[ particle_seeds ].z = (float)site_k;
		  particle_seed[ particle_seeds ].inlet_id =
		    (site_data & BOUNDARY_ID_MASK) >> BOUNDARY_ID_SHIFT;
		  ++particle_seeds;
		}
	}
  shared_vs = 0;
  
  for (m = 0; m < neigh_procs; m++)
    {
      shared_vs += neigh_proc[ m ].send_vs;
    }
  if (shared_vs > 0)
    {
      s_to_send = new short int[3 * shared_vs];
      s_to_recv = new short int[3 * shared_vs];
      
      v_to_send = new float[3 * shared_vs];
      v_to_recv = new float[3 * shared_vs];
    }
  shared_vs = 0;
  
  for (m = 0; m < neigh_procs; m++)
    {
      neigh_proc[ m ].s_to_send = &s_to_send[ shared_vs*3 ];
      neigh_proc[ m ].s_to_recv = &s_to_recv[ shared_vs*3 ];
      
      neigh_proc[ m ].v_to_send = &v_to_send[ shared_vs*3 ];
      neigh_proc[ m ].v_to_recv = &v_to_recv[ shared_vs*3 ];
      
      shared_vs += neigh_proc[ m ].send_vs;
    }
  for (m = 0; m < neigh_procs; m++)
    {
      neigh_proc[ m ].send_vs = 0;
    }
  particles_to_send_max = 1000;
  particles_to_recv_max = 1000;
  
  for (m = 0; m < neigh_procs; m++)
    {
      neigh_proc[ m ].p_to_send = new float[5 * particles_to_send_max];
      neigh_proc[ m ].p_to_recv = new float[5 * particles_to_recv_max];
    }
  
  req = new MPI_Request[2 * net->procs];
  
  from_proc_id_to_neigh_proc_index = new short int[net->procs];
  
  for (m = 0; m < net->procs; m++)
    {
      from_proc_id_to_neigh_proc_index[ m ] = -1;
    }
  for (m = 0; m < neigh_procs; m++)
    {
      from_proc_id_to_neigh_proc_index[ neigh_proc[m].id ] = m;
    }
  
  counter = 0;
  
  for (n = 0; n < blocks; n++)
    {
      if (velocity_field[ n ].vel_site_data == NULL) continue;
      
      for (m = 0; m < sites_in_a_block; m++)
	{
	  velocity_field[ n ].vel_site_data[ m ].counter = counter;
	}
      if (net->map_block[ n ].site_data == NULL) continue;
      
      for (m = 0; m < sites_in_a_block; m++)
	{
	  velocity_field[ n ].vel_site_data[ m ].site_id = net->map_block[ n ].site_data[ m ];
	}
    }
  procs = net->procs;
}


void streaklineDrawer::slRestart ()
{
  particles = 0;
}


void streaklineDrawer::slCommunicateSiteIds ()
{
#ifndef NOMPI
  int m;
  
  
  for (m = 0; m < neigh_procs; m++)
    {
      MPI_Irecv (&neigh_proc[ m ].recv_vs, 1, MPI_INT, neigh_proc[ m ].id,
		 30, MPI_COMM_WORLD, &req[ procs + neigh_proc[m].id ]);
      MPI_Isend (&neigh_proc[ m ].send_vs, 1, MPI_INT, neigh_proc[ m ].id,
		 30, MPI_COMM_WORLD, &req[ neigh_proc[m].id ]);
      MPI_Wait (&req[ neigh_proc[m].id ], status);
    }
  for (m = 0; m < neigh_procs; m++)
    {
      MPI_Wait (&req[ procs + neigh_proc[m].id ], status);
      
      if (neigh_proc[ m ].recv_vs > 0)
	{
	  MPI_Irecv (neigh_proc[ m ].s_to_recv, neigh_proc[ m ].recv_vs * 3, MPI_SHORT,
		     neigh_proc[ m ].id, 40, MPI_COMM_WORLD, &req[ procs + neigh_proc[m].id ]);
	}
      if (neigh_proc[ m ].send_vs > 0)
	{
	  MPI_Isend (neigh_proc[ m ].s_to_send, neigh_proc[ m ].send_vs * 3, MPI_SHORT,
		     neigh_proc[ m ].id, 40, MPI_COMM_WORLD, &req[ neigh_proc[m].id ]);
	  MPI_Wait (&req[ neigh_proc[m].id ], status);
	}
    }
  for (m = 0; m < neigh_procs; m++)
    {
      if (neigh_proc[ m ].recv_vs > 0)
	{
	  MPI_Wait (&req[ procs + neigh_proc[m].id ], status);
	}
    }
#endif // NOMPI
}


void streaklineDrawer::slCommunicateVelocities ()
{
#ifndef NOMPI
  int site_i, site_j, site_k;
  int neigh_i, neigh_j, neigh_k;
  int m, n;
  
  VelSiteData *vel_site_data_p;
  
  
  for (m = 0; m < neigh_procs; m++)
    {
      if (neigh_proc[ m ].send_vs > 0)
	{
	  MPI_Irecv (neigh_proc[ m ].v_to_recv, neigh_proc[ m ].send_vs * 3, MPI_FLOAT,
		     neigh_proc[ m ].id, 30, MPI_COMM_WORLD, &req[ procs + neigh_proc[m].id ]);
	}
    }
  for (m = 0; m < neigh_procs; m++)
    {
      for (n = 0; n < neigh_proc[ m ].recv_vs; n++)
	{
	  site_i = neigh_proc[ m ].s_to_recv[ 3*n+0 ];
	  site_j = neigh_proc[ m ].s_to_recv[ 3*n+1 ];
	  site_k = neigh_proc[ m ].s_to_recv[ 3*n+2 ];
	  
	  vel_site_data_p = slVelSiteDataPointer (site_i, site_j, site_k);
	  
	  if (vel_site_data_p != NULL)
	    {
	      neigh_proc[ m ].v_to_send[ 3*n+0 ] = vel_site_data_p->vx;
	      neigh_proc[ m ].v_to_send[ 3*n+1 ] = vel_site_data_p->vy;
	      neigh_proc[ m ].v_to_send[ 3*n+2 ] = vel_site_data_p->vz;
	    }
	  else
	    {
	      neigh_proc[ m ].v_to_send[ 3*n+0 ] = 0.;
	      neigh_proc[ m ].v_to_send[ 3*n+1 ] = 0.;
	      neigh_proc[ m ].v_to_send[ 3*n+2 ] = 0.;
	    }
	}
      if (neigh_proc[ m ].recv_vs > 0)
	{
	  MPI_Isend (neigh_proc[ m ].v_to_send, neigh_proc[ m ].recv_vs * 3, MPI_FLOAT,
		     neigh_proc[ m ].id, 30, MPI_COMM_WORLD, &req[ neigh_proc[m].id ]);
	  MPI_Wait (&req[ neigh_proc[m].id ], status);
	}
    }
  for (m = 0; m < neigh_procs; m++)
    {
      if (neigh_proc[ m ].send_vs <= 0) continue;
      
      MPI_Wait (&req[ procs + neigh_proc[m].id ], status);
      
      for (n = 0; n < neigh_proc[ m ].send_vs; n++)
	{
	  neigh_i = neigh_proc[ m ].s_to_send[ 3*n+0 ];
	  neigh_j = neigh_proc[ m ].s_to_send[ 3*n+1 ];
	  neigh_k = neigh_proc[ m ].s_to_send[ 3*n+2 ];
	  
	  vel_site_data_p = slVelSiteDataPointer (neigh_i, neigh_j, neigh_k);
	  
	  if (vel_site_data_p != NULL)
	    {
	      vel_site_data_p->vx = neigh_proc[ m ].v_to_recv[ 3*n+0 ];
	      vel_site_data_p->vy = neigh_proc[ m ].v_to_recv[ 3*n+1 ];
	      vel_site_data_p->vz = neigh_proc[ m ].v_to_recv[ 3*n+2 ];
	    }
	}
    }
  for (m = 0; m < neigh_procs; m++)
    {
      neigh_proc[ m ].send_vs = 0;
    }
#endif // NOMPI
}


void streaklineDrawer::slUpdateVelField (int stage_id, Net *net)
{
  float v[2][2][2][3], interp_v[3];
  float vel;
  
  int particles_temp;
  int is_interior;
  int n;
  
  
  particles_temp = particles;
  
  for (n = particles_temp - 1; n >= 0; n--)
    {
      slLocalVelField (n, v, &is_interior, net);
      
      if (stage_id == 0 && !is_interior) continue;
      
      slParticleVelocity (&particle[n], v, interp_v);
      vel = interp_v[0]*interp_v[0] + interp_v[1]*interp_v[1] + interp_v[2]*interp_v[2];
      
      if (vel > 1.0F)
	{
	  particle[n].vel = 1.0F;
	  particle[n].vx = interp_v[0] / sqrtf(vel);
	  particle[n].vy = interp_v[1] / sqrtf(vel);
	  particle[n].vz = interp_v[2] / sqrtf(vel);
	}
      else if (vel > 1.0e-8)
	{
	  particle[n].vel = sqrtf(vel);
	  particle[n].vx = interp_v[0];
	  particle[n].vy = interp_v[1];
	  particle[n].vz = interp_v[2];
	}
      else
      	{
      	  slDeleteParticle (n);
      	}
    }
}


void streaklineDrawer::slUpdateParticles ()
{
  for (int n = 0; n < particles; n++)
    {
      // particle coords updating (dt = 1)
      particle[n].x += particle[n].vx;
      particle[n].y += particle[n].vy;
      particle[n].z += particle[n].vz;
    }
}


void streaklineDrawer::slCommunicateParticles (Net *net)
{
#ifndef NOMPI
  int site_i, site_j, site_k;
  int m, n;
  int particles_temp;
  
  VelSiteData *vel_site_data_p;
  
  
  for (m = 0; m < neigh_procs; m++)
    {
      MPI_Irecv (&neigh_proc[ m ].recv_ps, 1, MPI_INT, neigh_proc[ m ].id,
		 30, MPI_COMM_WORLD, &req[ procs + neigh_proc[m].id ]);
    }
  for (m = 0; m < neigh_procs; m++)
    {
      neigh_proc[ m ].send_ps = 0;
    }
  particles_temp = particles;
  
  for (n = particles_temp - 1; n >= 0; n--)
    {
      site_i = (int)particle[ n ].x;
      site_j = (int)particle[ n ].y;
      site_k = (int)particle[ n ].z;
      
      vel_site_data_p = slVelSiteDataPointer (site_i, site_j, site_k);
      
      if (vel_site_data_p == NULL ||
	  vel_site_data_p->proc_id == net->id ||
	  vel_site_data_p->proc_id == -1)
	{
      	  continue;
      	}
      m = from_proc_id_to_neigh_proc_index[ vel_site_data_p->proc_id ];
      
      if (neigh_proc[ m ].send_ps == particles_to_send_max)
	{
	  particles_to_send_max *= 2;
	  neigh_proc[ m ].p_to_send =
	    (float *)realloc(neigh_proc[ m ].p_to_send, sizeof(float) * 5 * particles_to_send_max);
	}
      neigh_proc[ m ].p_to_send[ 5*neigh_proc[m].send_ps+0 ] = particle[ n ].x;
      neigh_proc[ m ].p_to_send[ 5*neigh_proc[m].send_ps+1 ] = particle[ n ].y;
      neigh_proc[ m ].p_to_send[ 5*neigh_proc[m].send_ps+2 ] = particle[ n ].z;
      neigh_proc[ m ].p_to_send[ 5*neigh_proc[m].send_ps+3 ] = particle[ n ].vel;
      neigh_proc[ m ].p_to_send[ 5*neigh_proc[m].send_ps+4 ] = particle[ n ].inlet_id + 0.1;
      ++neigh_proc[ m ].send_ps;
      
      slDeleteParticle (n);
    }
  for (m = 0; m < neigh_procs; m++)
    {
      MPI_Isend (&neigh_proc[ m ].send_ps, 1, MPI_INT, neigh_proc[ m ].id,
		 30, MPI_COMM_WORLD, &req[ neigh_proc[m].id ]);
      MPI_Wait (&req[ neigh_proc[m].id ], net->status);
    }
  for (m = 0; m < neigh_procs; m++)
    {
      MPI_Wait (&req[ procs + neigh_proc[m].id ], net->status);
    }
  
  for (m = 0; m < neigh_procs; m++)
    {
      if (neigh_proc[ m ].send_ps > 0)
	{
	  MPI_Isend (neigh_proc[ m ].p_to_send, neigh_proc[ m ].send_ps * 5, MPI_FLOAT,
		     neigh_proc[ m ].id, 40, MPI_COMM_WORLD, &req[ neigh_proc[m].id ]);
	  MPI_Wait (&req[ neigh_proc[m].id ], net->status);
	}
    }
  for (m = 0; m < neigh_procs; m++)
    {
      if (neigh_proc[ m ].recv_ps > 0)
	{
	  if (neigh_proc[ m ].recv_ps > particles_to_recv_max)
	    {
	      particles_to_recv_max *= 2;
	      particles_to_recv_max = UtilityFunctions::max(particles_to_recv_max, neigh_proc[ m ].recv_ps);
	      neigh_proc[ m ].p_to_recv =
	      	(float *)realloc(neigh_proc[ m ].p_to_recv, sizeof(float) * 5 * particles_to_recv_max);
	    }
	  MPI_Irecv (neigh_proc[ m ].p_to_recv, neigh_proc[ m ].recv_ps * 5, MPI_FLOAT,
		     neigh_proc[ m ].id, 40, MPI_COMM_WORLD, &req[ procs + neigh_proc[m].id ]);
	  MPI_Wait (&req[ procs + neigh_proc[m].id ], net->status);
	  
	  for (n = 0; n < neigh_proc[ m ].recv_ps; n++)
	    {
	      slCreateParticle (neigh_proc[ m ].p_to_recv[ 5*n+0 ],
				neigh_proc[ m ].p_to_recv[ 5*n+1 ],
				neigh_proc[ m ].p_to_recv[ 5*n+2 ],
				neigh_proc[ m ].p_to_recv[ 5*n+3 ],
				(int)neigh_proc[ m ].p_to_recv[ 5*n+4 ]);
	    }
	}
    }
#endif // NOMPI
}


void streaklineDrawer::render ()
{
  float screen_max[2];
  float scale[2];
  float p1[3], p2[3];
  
  int pixels_x, pixels_y;
  int i, j;
  int n;
  
  ColPixel col_pixel;
  
  
  pixels_x = screen.pixels_x;
  pixels_y = screen.pixels_y;
  
  screen_max[0] = screen.max_x;
  screen_max[1] = screen.max_y;
  
  scale[0] = screen.scale_x;
  scale[1] = screen.scale_y;
  
  for (n = 0; n < particles; n++)
    {
      p1[0] = particle[n].x - (float)(sites_x>>1);
      p1[1] = particle[n].y - (float)(sites_y>>1);
      p1[2] = particle[n].z - (float)(sites_z>>1);
      
      visProject (p1, p2);
      
      p2[0] = (int)(scale[0] * (p2[0] + screen_max[0]));
      p2[1] = (int)(scale[1] * (p2[1] + screen_max[1]));
      
      i = (int)p2[0];
      j = (int)p2[1];
      
      if (!(i < 0 || i >= pixels_x ||
	    j < 0 || j >= pixels_y))
	{
	  col_pixel.particle_vel      = particle[n].vel;
	  col_pixel.particle_z        = p2[2];
	  col_pixel.particle_inlet_id = particle[n].inlet_id;
	  col_pixel.i                 = PixelId(i, j) | STREAKLINE;
	  
	  visWritePixel (&col_pixel);
	}
    }
}


void streaklineDrawer::slStreakLines (int time_steps, int time_steps_per_cycle, Net *net)
{
  if (!is_bench)
    {
      int particle_creation_period = UtilityFunctions::max(1, (int)(time_steps_per_cycle / 5000.0F));
      
      if (time_steps % (int)(time_steps_per_cycle / vis_streaklines_per_pulsatile_period) <=
	  (vis_streakline_length/100.0F) * (time_steps_per_cycle / vis_streaklines_per_pulsatile_period) &&
	  time_steps % particle_creation_period == 0)
	{
	  slCreateSeedParticles ();
	}
    }
  else
    {
      if (time_steps % 10 == 0)
	{
	  slCreateSeedParticles ();
	}
    }
  ++counter;
  
  slUpdateVelField (0, net);
  slCommunicateSiteIds ();
  slCommunicateVelocities ();
  slUpdateVelField (1, net);
  slUpdateParticles ();
  slCommunicateParticles (net);
}


void streaklineDrawer::slEnd ()
{
  int m;
  
  
  free(from_proc_id_to_neigh_proc_index);
  
  free(req);
  
  for (m = 0; m < neigh_procs; m++)
    {
      free(neigh_proc[ m ].p_to_recv);
      free(neigh_proc[ m ].p_to_send);
    }
  
  if (shared_vs > 0)
    {
      free(v_to_recv);
      free(v_to_send);
      
      free(s_to_recv);
      free(s_to_send);
    }
  for (m = 0; m < blocks; m++)
    {
      if (velocity_field[ m ].vel_site_data != NULL)
	{
	  free(velocity_field[ m ].vel_site_data);
	}
    }
  free(velocity_field);
  
  free(particle_seed);
  
  free(particle);
}