// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

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
                           const VisSettings& iVisSettings,
                           const lb::MacroscopicPropertyCache& propertyCache);
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
          const lb::MacroscopicPropertyCache& propertyCache;

          std::map<proc_t, NeighbouringProcessor> neighbouringProcessors;
          ParticleManager particleManager;
          VelocityField velocityField;

          std::vector<Particle> particleSeeds;
      };
    }
  }
}

#endif // HEMELB_VIS_STREAKLINEDRAWER_STREAKLINEDRAWER_H
