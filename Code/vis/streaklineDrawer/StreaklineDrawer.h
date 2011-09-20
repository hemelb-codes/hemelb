#ifndef HEMELB_VIS_STREAKLINEDRAWER_H
#define HEMELB_VIS_STREAKLINEDRAWER_H

#include <vector>

#include "constants.h"
#include "mpiInclude.h"

#include "geometry/LatticeData.h"
#include "topology/NetworkTopology.h"

#include "vis/Screen.h"
#include "vis/Viewpoint.h"
#include "vis/VisSettings.h"

namespace hemelb
{
  namespace vis
  {
    namespace streaklinedrawer
    {

      //TODO Some of this class's members could be combined with the topology-aware classes
      // (and obv be initialised when they are).

      /**
       * Class that controls the drawing of streaklines - lines that trace
       * the path of an imaginary particle were it dropped into the fluid.
       */
      class StreaklineDrawer
      {
      public:
        // Constructor and destructor.
        StreaklineDrawer(const geometry::LatticeData& iLatDat,
                         Screen& iScreen,
                         const Viewpoint& iViewpoint,
                         const VisSettings& iVisSettings);
        ~StreaklineDrawer();

        // Method to reset streakline drawer
        void Restart();

        // Drawing methods.
        void StreakLines(unsigned long time_steps,
                         unsigned long time_steps_per_cycle,
                         const geometry::LatticeData& iLatDat);
	
        void render(const geometry::LatticeData& iLatDat);

      private:

        // Struct for a particle dropped into the fluid.
        struct Particle
        {
	  float x, y, z;
	  float vx, vy, vz;
	  float vel;
	  unsigned int inlet_id;
        };

        // Struct for information about the velocity field at some point.
        struct VelSiteData
        {
	public:
	  VelSiteData()
	    {
	      proc_id = -1;
	      counter = 0;
	    }

	  proc_t proc_id;
	  site_t counter, site_id;
	  float vx, vy, vz;
        };

	// Struct for one processor to hold information about its neighbouring processors.
        typedef std::vector<float> FloatVector;
        struct NeighProc
        {
	  proc_t id;
	  site_t send_ps, recv_ps;
	  site_t send_vs, recv_vs;

	  FloatVector p_to_send, p_to_recv;
	  float *v_to_send, *v_to_recv;

	  site_t *s_to_send, *s_to_recv;
        };

        // Counter keeps track of the number of VelSiteDatas created
        site_t counter;

        // Variables for tracking the actual numbers of particles, and the maximum number
        //(i.e. the number for which memory has been allocated).
        unsigned int nParticles;
        unsigned int nParticleSeeds;
        unsigned int particles_to_send_max, particles_to_recv_max;

        // Variables for counting the processors involved etc.
        site_t shared_vs;
        proc_t procs;

	//Vector containing VelocityFields
	std::vector<std::vector<VelSiteData> > velocity_field;
	
        std::vector<Particle> mParticleVec;
        std::vector<Particle> mParticleSeedVec;

        // Arrays for communicating between processors.
        float *v_to_send, *v_to_recv;
        site_t *s_to_send, *s_to_recv;
        proc_t *from_proc_id_to_neigh_proc_index;

        std::vector<NeighProc> mNeighProcs;

        // Require these for inter-processor comms.
        MPI_Request *req;

        // Private functions for the creation / deletion of particles.
        void createSeedParticles();
        void createParticle(float x, float y, float z, float vel, int inlet_id);
        void deleteParticle(unsigned int p_index);

        // Private functions for initialising the velocity field.
        void initializeVelFieldBlock(const geometry::LatticeData& iLatDat,
                                     site_t site_i,
                                     site_t site_j,
                                     site_t site_k,
                                     proc_t proc_id);
	
        VelSiteData *velSiteDataPointer(const geometry::LatticeData& iLatDat,
                                        site_t site_i,
                                        site_t site_j,
                                        site_t site_k);
        void particleVelocity(Particle *particle_p, float v[2][2][2][3], float interp_v[3]);
       
	void localVelField(int p_index,
                           float v[2][2][2][3],
                           int *is_interior,
                           const geometry::LatticeData& iLatDat);

        // Private functions for updating the velocity field and the particles in it.
        void updateVelField(int stage_id, const geometry::LatticeData& iLatDat);
        void updateParticles();

        // Private functions for inter-proc communication.
        void communicateSiteIds();
        void communicateVelocities(const geometry::LatticeData& iLatDat);
        void communicateParticles(const geometry::LatticeData& iLatDat);

	Screen& mScreen;
        const Viewpoint& mViewpoint;
        const VisSettings& mVisSettings;
      };
    }
  }
}

#endif // HEMELB_VIS_STREAKLINEDRAWER_H
