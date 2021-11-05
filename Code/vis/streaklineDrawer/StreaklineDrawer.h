// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_VIS_STREAKLINEDRAWER_STREAKLINEDRAWER_H
#define HEMELB_VIS_STREAKLINEDRAWER_STREAKLINEDRAWER_H

#include <vector>
#include <map>

#include "constants.h"
#include "net/mpi.h"

#include "geometry/LatticeData.h"
#include "net/IOCommunicator.h"
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
          StreaklineDrawer(const geometry::LatticeData& iLatDat, const Screen& iScreen,
                           const Viewpoint& iViewpoint, const VisSettings& iVisSettings,
                           const lb::MacroscopicPropertyCache& propertyCache,
                           const net::MpiCommunicator& comms);
          ~StreaklineDrawer();

          // Method to reset streakline drawer
          void Restart();

          // Drawing methods.
          void ProgressStreaklines(unsigned long time_steps, unsigned long total_time_steps);
          PixelSet<StreakPixel>* Render();

        private:
          // Function for updating the velocity field and the particles in it.
          void UpdateVelocityFieldForAllParticlesAndPrune();
          void UpdateVelocityFieldForCommunicatedSites();

          // Private functions for the creation / deletion of particles.
          void ChooseSeedParticles();
          void CreateParticlesFromSeeds();

          // Private functions for inter-proc communication.
          void CommunicateSiteIds();
          void CommunicateVelocities();
          void WorkOutVelocityDataNeededForParticles();

          const geometry::LatticeData& latDat;
          const Screen& screen;
          const Viewpoint& viewpoint;
          const VisSettings& visSettings;
          //const lb::MacroscopicPropertyCache& propertyCache;

          std::map<proc_t, NeighbouringProcessor> neighbouringProcessors;
          ParticleManager particleManager;
          VelocityField velocityField;

          std::vector<Particle> particleSeeds;
          net::Net* streakNet;
      };
    }
  }
}

#endif // HEMELB_VIS_STREAKLINEDRAWER_STREAKLINEDRAWER_H
