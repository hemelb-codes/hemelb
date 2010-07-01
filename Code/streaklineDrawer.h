#ifndef __streaklineDrawer_h_
#define __streaklineDrawer_h_

#include "constants.h"
#include "mpiInclude.h"
#include "visualisationControl.h"
#include "net.h"

class streaklineDrawer : public visualisationLayer
{
  private:

    struct Particle
    {
      float x, y, z;
      float vx, vy, vz;
      float vel;
      int inlet_id;
    };

    struct VelSiteData
    {
      int proc_id, counter, site_id;
      float vx, vy, vz;
    };

    struct VelocityField
    {
      VelSiteData *vel_site_data;
    };

    struct NeighProc
    {
      int id;
      int send_ps, recv_ps;
      int send_vs, recv_vs;
      
      float *p_to_send, *p_to_recv;
      float *v_to_send, *v_to_recv;
    
      short int *s_to_send, *s_to_recv;
    };

    int counter;
    int particles, particles_max;
    int particle_seeds, particle_seeds_max;
    int particles_to_send_max, particles_to_recv_max;
    int neigh_procs;
    int shared_vs;
    int procs;
  
    VelocityField *velocity_field;
    Particle *particle;
    Particle *particle_seed;
  
    float *v_to_send, *v_to_recv;
    short int *s_to_send, *s_to_recv;
    short int *from_proc_id_to_neigh_proc_index;
  
    NeighProc neigh_proc[NEIGHBOUR_PROCS_MAX];
  
    #ifndef NOMPI
      MPI_Status status[4];
      MPI_Request *req;
    #endif

    void slCreateSeedParticles ();
    void slCreateParticle (float x, float y, float z, float vel, int inlet_id);
    void slDeleteParticle (int p_index);
  
    void slInitializeVelFieldBlock(int site_i, int site_j, int site_k, int proc_id);
    VelSiteData *slVelSiteDataPointer (int site_i, int site_j, int site_k);
    void slParticleVelocity (Particle *particle_p, float v[2][2][2][3], float interp_v[3]);
    void slLocalVelField (int p_index, float v[2][2][2][3], int *is_interior, Net *net);
  
    void slUpdateVelField (int stage_id, Net *net);
    void slUpdateParticles ();

    void slCommunicateSiteIds ();
    void slCommunicateVelocities ();
    void slCommunicateParticles (Net *net);


  public:
    void slStreakLines (int time_steps, int time_steps_per_cycle, Net *net);
    virtual void render ();
    void slInit (Net *net);
    void slRestart ();
    void slEnd ();
};

#endif //__streaklineDrawer_h_