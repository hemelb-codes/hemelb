#ifndef HEMELB_VIS_STREAKLINEDRAWER_H
#define HEMELB_VIS_STREAKLINEDRAWER_H

#include <vector>

#include "constants.h"
#include "mpiInclude.h"

#include "geometry/LatticeData.h"
#include "topology/NetworkTopology.h"

#include "vis/streaklineDrawer/NeighProc.h"
#include "vis/streaklineDrawer/Particle.h"
#include "vis/streaklineDrawer/Particles.h"
#include "vis/streaklineDrawer/VelocityField.h"
#include "vis/streaklineDrawer/VelocitySiteData.h"
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
	// Private functions for updating the velocity field and the particles in it.
        void updateVelField(int stage_id, const geometry::LatticeData& iLatDat);
        void updateParticles();
	
     public:
	// Variables for counting the processors involved etc.
        site_t shared_vs;
        proc_t procs;


      
        std::vector<Particle> mParticleSeeds;

      private:



        // Arrays for communicating between processors.
        float *v_to_send, *v_to_recv;
        site_t *s_to_send, *s_to_recv;
	
      public:
        proc_t *from_proc_id_to_neigh_proc_index;


        std::vector<NeighProc> mNeighProcs;

      private:      
        // Require these for inter-processor comms.
        MPI_Request *req;

        // Private functions for the creation / deletion of particles.
        void createSeedParticles();
    

        // Private functions for inter-proc communication.
        void communicateSiteIds();
        void communicateVelocities(const geometry::LatticeData& iLatDat);
       
	VelocityField mVelocityField;

	Particles mParticles;

	Screen& mScreen;
        const Viewpoint& mViewpoint;
        const VisSettings& mVisSettings;
      };
    }
  }
}

#endif // HEMELB_VIS_STREAKLINEDRAWER_H
