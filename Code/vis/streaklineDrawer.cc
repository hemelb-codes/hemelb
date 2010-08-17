// TODO: Mixture of malloc and new is a bad smell, and just using new / delete would be much nicer. 
// Vectors might be a good way of getting round the need to use realloc.
#include <stdlib.h>
#include <math.h>

#include "vis/StreaklineDrawer.h"
// TODO this could probably be reduced to the net class and some visualisation class.
#include "vis/rt.h"
#include "utilityFunctions.h"

using namespace vis;

// Function to initialise the velocity field at given coordinates.
void StreaklineDrawer::initializeVelFieldBlock (int site_i, int site_j, int site_k, int proc_id)
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

// Returns the velocity site data for a given index, or NULL if the index isn't valid / has
// no data.
StreaklineDrawer::VelSiteData *StreaklineDrawer::velSiteDataPointer (int site_i, int site_j, int site_k)
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

// Interpolates a velocity field to get the velocity at the position of a particle.
void StreaklineDrawer::particleVelocity (Particle *particle_p, float v[2][2][2][3], float interp_v[3])
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

// Create a particle with given position, velocity and inlet_id.
void StreaklineDrawer::createParticle (float x, float y, float z, float vel, int inlet_id)
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

// Delete the particle at given index. Do something a bit budget to ensure that 
// the particles remain in the first <particles> elements of an array,
void StreaklineDrawer::deleteParticle (int p_index)
{
  if (--particles <= 0) return;
  
  // its data are replaced with those of the last particle;
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

// Create seed particles to begin the streaklines.
void StreaklineDrawer::createSeedParticles ()
{
  for (int n = 0; n < particle_seeds; n++)
    {
      createParticle (particle_seed[n].x,
			particle_seed[n].y,
			particle_seed[n].z,
			0.0F,
			particle_seed[n].inlet_id);
    }
}

// Populate the matrix v with all the velocity field data at each index.
void StreaklineDrawer::localVelField (int p_index, float v[2][2][2][3], int *is_interior, Net *net)
{
  double vx, vy, vz;
  
  VelSiteData *vel_site_data_p;
    
  int site_i = (int)particle[ p_index ].x;
  int site_j = (int)particle[ p_index ].y;
  int site_k = (int)particle[ p_index ].z;
  
  *is_interior = 1;
  
  int c1Plusc2 = check_conv ? 45 : 15;

  for (int i = 0; i < 2; i++)
    {
      int neigh_i = site_i + i;
      
      for (int j = 0; j < 2; j++)
	{
	  int neigh_j = site_j + j;
	  
	  for (int k = 0; k < 2; k++)
	    {
	      int neigh_k = site_k + k;
	      
	      vel_site_data_p = velSiteDataPointer (neigh_i, neigh_j, neigh_k);
	      
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
		  // This means that the local velocity has already been
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
		  double density;

		  lbmDensityAndVelocity (&f_old[ vel_site_data_p->site_id*c1Plusc2 ],
					 &density, &vx, &vy, &vz);
		  
		  v[i][j][k][0] = vel_site_data_p->vx = vx / density;
		  v[i][j][k][1] = vel_site_data_p->vy = vy / density;
		  v[i][j][k][2] = vel_site_data_p->vz = vz / density;
		}
	      else
		{
		  vel_site_data_p->counter = counter;
		  
		  int m = from_proc_id_to_neigh_proc_index[ vel_site_data_p->proc_id ];
		  
		  neigh_proc[m].s_to_send[ 3*neigh_proc[m].send_vs+0 ] = neigh_i;
		  neigh_proc[m].s_to_send[ 3*neigh_proc[m].send_vs+1 ] = neigh_j;
		  neigh_proc[m].s_to_send[ 3*neigh_proc[m].send_vs+2 ] = neigh_k;
		  ++neigh_proc[m].send_vs;
		}
	    }
	}
    }
}

// Constructor, populating fields from a Net object.
StreaklineDrawer::StreaklineDrawer (Net *net)
{
  int m, mm, n;

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
  
  for (int i = 0; i < sites_x; i += block_size)
    for (int j = 0; j < sites_y; j += block_size)
      for (int k = 0; k < sites_z; k += block_size)
	{
	  map_block_p = &net->map_block[ ++n ];
	  
	  if (map_block_p->site_data == NULL) continue;
	  
	  proc_block_p = &net->proc_block[ n ];
	  
	  m = -1;
	  
	  for (int site_i = i; site_i < i + block_size; site_i++)
	    for (int site_j = j; site_j < j + block_size; site_j++)
	      for (int site_k = k; site_k < k + block_size; site_k++)
		{
		  if (proc_block_p->proc_id[ ++m ] != net->id) continue;
		  
		  for (int neigh_i = util::max(0, site_i-1); neigh_i <= util::min(sites_x-1, site_i+1); neigh_i++)
		    for (int neigh_j = util::max(0, site_j-1); neigh_j <= util::min(sites_y-1, site_j+1); neigh_j++)
		      for (int neigh_k = util::max(0, site_k-1); neigh_k <= util::min(sites_z-1, site_k+1); neigh_k++)
			{
			  neigh_proc_id = net->netProcIdPointer (neigh_i, neigh_j, neigh_k);
			  
			  if (neigh_proc_id == NULL || *neigh_proc_id == (1 << 30))
			    {
			      continue;
			    }
			  initializeVelFieldBlock (neigh_i, neigh_j, neigh_k, *neigh_proc_id);
			  
			  if (*neigh_proc_id == net->id) continue;
			  
			  vel_site_data_p = velSiteDataPointer (neigh_i, neigh_j, neigh_k);
			  
			  if (vel_site_data_p->counter == counter) continue;
			  
			  vel_site_data_p->counter = counter;
			  
                          bool seenSelf = false;
			  for (mm = 0; mm < neigh_procs && !seenSelf; mm++)
			    {
			      if (*neigh_proc_id == neigh_proc[ mm ].id)
				{
				  seenSelf = true;
				  ++neigh_proc[ mm ].send_vs;
				}
			    }
			  if (seenSelf) continue;
			  
			  if (neigh_procs == NEIGHBOUR_PROCS_MAX)
			    {
			      printf (" too many inter processor neighbours in streakline constructor()\n");
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

// Reset the streakline drawer.
void StreaklineDrawer::restart ()
{
  particles = 0;
}

// Communicate site ids to other processors.
void StreaklineDrawer::communicateSiteIds ()
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

// Communicate velocities to other processors.
void StreaklineDrawer::communicateVelocities ()
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
	  
	  vel_site_data_p = velSiteDataPointer (site_i, site_j, site_k);
	  
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
	  
	  vel_site_data_p = velSiteDataPointer (neigh_i, neigh_j, neigh_k);
	  
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

// Update the velocity field.
void StreaklineDrawer::updateVelField (int stage_id, Net *net)
{
  float v[2][2][2][3], interp_v[3];
  float vel;
  
  int particles_temp;
  int is_interior;
  int n;
  
  
  particles_temp = particles;
  
  for (n = particles_temp - 1; n >= 0; n--)
    {
      localVelField (n, v, &is_interior, net);
      
      if (stage_id == 0 && !is_interior) continue;
      
      particleVelocity (&particle[n], v, interp_v);
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
      	  deleteParticle (n);
      	}
    }
}

// Update the particles.
void StreaklineDrawer::updateParticles ()
{
  for (int n = 0; n < particles; n++)
    {
      // particle coords updating (dt = 1)
      particle[n].x += particle[n].vx;
      particle[n].y += particle[n].vy;
      particle[n].z += particle[n].vz;
    }
}

// Communicate that particles current state to other processors.
void StreaklineDrawer::communicateParticles (Net *net)
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
      
      vel_site_data_p = velSiteDataPointer (site_i, site_j, site_k);
      
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
      
      deleteParticle (n);
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
	      particles_to_recv_max = util::max(particles_to_recv_max, neigh_proc[ m ].recv_ps);
	      neigh_proc[ m ].p_to_recv =
	      	(float *)realloc(neigh_proc[ m ].p_to_recv, sizeof(float) * 5 * particles_to_recv_max);
	    }
	  MPI_Irecv (neigh_proc[ m ].p_to_recv, neigh_proc[ m ].recv_ps * 5, MPI_FLOAT,
		     neigh_proc[ m ].id, 40, MPI_COMM_WORLD, &req[ procs + neigh_proc[m].id ]);
	  MPI_Wait (&req[ procs + neigh_proc[m].id ], net->status);
	  
	  for (n = 0; n < neigh_proc[ m ].recv_ps; n++)
	    {
	      createParticle (neigh_proc[ m ].p_to_recv[ 5*n+0 ],
				neigh_proc[ m ].p_to_recv[ 5*n+1 ],
				neigh_proc[ m ].p_to_recv[ 5*n+2 ],
				neigh_proc[ m ].p_to_recv[ 5*n+3 ],
				(int)neigh_proc[ m ].p_to_recv[ 5*n+4 ]);
	    }
	}
    }
#endif // NOMPI
}

// Render the streaklines
void StreaklineDrawer::render ()
{
  float screen_max[2];
  float scale[2];
  float p1[3], p2[3];
  
  ColPixel col_pixel;
  
  int pixels_x = screen.pixels_x;
  int pixels_y = screen.pixels_y;
  
  screen_max[0] = screen.max_x;
  screen_max[1] = screen.max_y;
  
  scale[0] = screen.scale_x;
  scale[1] = screen.scale_y;
  
  for (int n = 0; n < particles; n++)
    {
      p1[0] = particle[n].x - (float)(sites_x>>1);
      p1[1] = particle[n].y - (float)(sites_y>>1);
      p1[2] = particle[n].z - (float)(sites_z>>1);
      
      visProject (p1, p2);
      
      p2[0] = (int)(scale[0] * (p2[0] + screen_max[0]));
      p2[1] = (int)(scale[1] * (p2[1] + screen_max[1]));
      
      int i = (int)p2[0];
      int j = (int)p2[1];
      
      if (!(i < 0 || i >= pixels_x ||
	    j < 0 || j >= pixels_y))
	{
	  col_pixel.particle_vel      = particle[n].vel;
	  col_pixel.particle_z        = p2[2];
	  col_pixel.particle_inlet_id = particle[n].inlet_id;
	  col_pixel.i                 = PixelId(i, j);
	  col_pixel.i.isStreakline    = true;
	  
	  visWritePixel (&col_pixel);
	}
    }
}

// Draw streaklines
void StreaklineDrawer::streakLines (int time_steps, int time_steps_per_cycle, Net *net)
{
  if (!is_bench)
    {
      int particle_creation_period = util::max(1, (int)(time_steps_per_cycle / 5000.0F));
      
      if (time_steps % (int)(time_steps_per_cycle / streaklines_per_pulsatile_period) <=
	  (streakline_length/100.0F) * (time_steps_per_cycle / streaklines_per_pulsatile_period) &&
	  time_steps % particle_creation_period == 0)
	{
	  createSeedParticles ();
	}
    }
  else
    {
      if (time_steps % 10 == 0)
	{
	  createSeedParticles ();
	}
    }
  ++counter;
  
  updateVelField (0, net);
  communicateSiteIds ();
  communicateVelocities ();
  updateVelField (1, net);
  updateParticles ();
  communicateParticles (net);
}

// Destructor
StreaklineDrawer::~StreaklineDrawer ()
{
  free(from_proc_id_to_neigh_proc_index);
  free(req);
  
  for (int m = 0; m < neigh_procs; m++)
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
  for (int m = 0; m < blocks; m++)
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
