#ifndef __vis_streaklineDrawer_h_
#define __vis_streaklineDrawer_h_

#include "constants.h"
#include "mpiInclude.h"
#include "vis/visualisationControl.h"
#include "net.h"

namespace vis {

  // Class that controls the drawing of streaklines - lines that trace the path of an
  // imaginary particle were it dropped into the fluid.
  class streaklineDrawer : public visualisationLayer
  {
  private:

    // Struct for a particle dropped into the fluid.
    // TODO: Get rid of inlet_id because I think it's never meaningfully used.
    struct Particle
    {
      float x, y, z;
      float vx, vy, vz;
      float vel;
      int inlet_id;
    };

    // Struct for information about the velocity field at some point.
    struct VelSiteData
    {
      int proc_id, counter, site_id;
      float vx, vy, vz;
    };

    // Struct to contain the whole velocity field.
    struct VelocityField
    {
      VelSiteData *vel_site_data;
    };

    // Struct for one processor to hold information about its neighbouring processors.
    struct NeighProc
    {
      int id;
      int send_ps, recv_ps;
      int send_vs, recv_vs;
      
      float *p_to_send, *p_to_recv;
      float *v_to_send, *v_to_recv;
    
      short int *s_to_send, *s_to_recv;
    };

    // Counter keeps track of the number of VelSiteDatas created
    int counter;

    // Variables for tracking the actual numbers of particles, and the maximum number 
    //(i.e. the number for which memory has been allocated).
    int particles, particles_max;
    int particle_seeds, particle_seeds_max;
    int particles_to_send_max, particles_to_recv_max;

    // Variables for counting the processors involved etc.
    int neigh_procs;
    int shared_vs;
    int procs;
  
    // Pointers to the structs.
    VelocityField *velocity_field;
    Particle *particle;
    Particle *particle_seed;
  
    // Arrays for communicating between processors.
    float *v_to_send, *v_to_recv;
    short int *s_to_send, *s_to_recv;
    short int *from_proc_id_to_neigh_proc_index;
  
    NeighProc neigh_proc[NEIGHBOUR_PROCS_MAX];
  
    // If using MPI, require these for inter-processor comms.
#ifndef NOMPI
    MPI_Status status[4];
    MPI_Request *req;
#endif

    // Private functions for the creation / deletion of particles.
    void createSeedParticles ();
    void createParticle (float x, float y, float z, float vel, int inlet_id);
    void deleteParticle (int p_index);

    // Private functions for initialising the velocity field.  
    void initializeVelFieldBlock(int site_i, int site_j, int site_k, int proc_id);
    VelSiteData *velSiteDataPointer (int site_i, int site_j, int site_k);
    void particleVelocity (Particle *particle_p, float v[2][2][2][3], float interp_v[3]);
    void localVelField (int p_index, float v[2][2][2][3], int *is_interior, Net *net);
  
    // Private functions for updating the velocity field and the particles in it.
    void updateVelField (int stage_id, Net *net);
    void updateParticles ();

    // Private functions for inter-proc communication.
    void communicateSiteIds ();
    void communicateVelocities ();
    void communicateParticles (Net *net);

  public:
    // Constructor and destructor.
    streaklineDrawer (Net *net);
    ~streaklineDrawer ();

    // Method to reset streakline drawer
    void restart ();

    // Drawing methods.
    void streakLines (int time_steps, int time_steps_per_cycle, Net *net);
    virtual void render ();
  };

}

#endif //__vis_streaklineDrawer_h_
