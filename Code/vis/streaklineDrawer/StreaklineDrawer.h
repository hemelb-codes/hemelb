#ifndef HEMELB_VIS_STREAKLINEDRAWER_STREAKLINEDRAWER_H
#define HEMELB_VIS_STREAKLINEDRAWER_STREAKLINEDRAWER_H

#include <vector>
#include <map>

#include "constants.h"
#include "mpiInclude.h"

#include "debug/Debugger.h"
#include "geometry/LatticeData.h"
#include "topology/NetworkTopology.h"
#include "vis/PixelSet.h"
#include "vis/PixelSetStore.h"
#include "vis/streaklineDrawer/NeighbouringProcessor.h"
#include "vis/streaklineDrawer/ParticleManager.h"
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
                           const Screen& iScreen,
                           const Viewpoint& iViewpoint,
                           const VisSettings& iVisSettings);
          ~StreaklineDrawer();

          // Method to reset streakline drawer
          void Restart();

          // Drawing methods.
          void StreakLines(unsigned long time_steps, unsigned long time_steps_per_cycle);

          PixelSet<StreakPixel>* Render();

          std::vector<Particle> particleSeeds;

        private:
          // Function for updating the velocity field and the particles in it.
          void UpdateVelocityFieldForAllParticles();
          void UpdateVelocityFieldForCommunicatedSites();

          void updateParticles();

          // Variables for counting the processors involved etc.
          proc_t procs;

          // Private functions for the creation / deletion of particles.
          void createSeedParticles();

          // Private functions for inter-proc communication.
          void CommunicateSiteIds();
          void CommunicateVelocities();

          std::map<proc_t, NeighbouringProcessor> neighbouringProcessors;

          const geometry::LatticeData& latDat;
          ParticleManager particleManager;
          const Screen& mScreen;
          VelocityField velocityField;
          const Viewpoint& mViewpoint;
          const VisSettings& mVisSettings;
      };
    }
  }
}

#endif // HEMELB_VIS_STREAKLINEDRAWER_STREAKLINEDRAWER_H
