#ifndef HEMELB_VIS_STREAKLINEDRAWER_H
#define HEMELB_VIS_STREAKLINEDRAWER_H

#include <vector>
#include <map>

#include "constants.h"
#include "mpiInclude.h"

#include "geometry/LatticeData.h"
#include "topology/NetworkTopology.h"
#include "vis/PixelSet.h"
#include "vis/PixelSetStore.h"
#include "vis/streaklineDrawer/NeighbouringProcessor.h"
#include "vis/streaklineDrawer/Particle.h"
#include "vis/streaklineDrawer/Particles.h"
#include "vis/streaklineDrawer/VelocityField.h"
#include "vis/streaklineDrawer/VelocitySiteData.h"
#include "vis/Screen.h"
#include "vis/streaklineDrawer/StreakPixel.h"
#include "vis/Viewpoint.h"
#include "vis/VisSettings.h"

namespace hemelb
{
  namespace vis
  {
    namespace streaklinedrawer
    {
      /**
       * Class that controls the drawing of streaklines - lines that trace
       * the path of an imaginary particle were it dropped into the fluid.
       */
      class StreaklineDrawer : public PixelSetStore<PixelSet<StreakPixel> >
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
          void StreakLines(unsigned long time_steps, unsigned long time_steps_per_cycle);

          PixelSet<StreakPixel>* Render();

          std::vector<Particle> mParticleSeeds;
          proc_t *from_proc_id_to_neigh_proc_index;

        private:
          // Functions for updating the velocity field and the particles in it.
          void updateVelField(int stage_id);
          void updateParticles();

          // Arrays for communicating between processors.
          float *v_to_send, *v_to_recv;
          site_t *s_to_send, *s_to_recv;

          // Variables for counting the processors involved etc.
          site_t shared_vs;
          proc_t procs;

          // Require these for inter-processor comms.
          MPI_Request *req;

          // Private functions for the creation / deletion of particles.
          void createSeedParticles();

          // Private functions for inter-proc communication.
          void communicateSiteIds();
          void communicateVelocities();

          std::vector<NeighbouringProcessor> mNeighbouringProcessors;
          const geometry::LatticeData& latDat;
          Particles mParticles;
          Screen& mScreen;
          VelocityField mVelocityField;
          const Viewpoint& mViewpoint;
          const VisSettings& mVisSettings;
      };
    }
  }
}

#endif // HEMELB_VIS_STREAKLINEDRAWER_H
