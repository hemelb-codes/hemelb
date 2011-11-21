#ifndef HEMELB_VIS_STREAKLINEDRAWER_PARTICLEMANAGER_H
#define HEMELB_VIS_STREAKLINEDRAWER_PARTICLEMANAGER_H

#include <vector>
#include <map>

#include "constants.h"
#include "mpiInclude.h"

#include "geometry/LatticeData.h"
#include "topology/NetworkTopology.h"

#include "vis/streaklineDrawer/NeighbouringProcessor.h"
#include "vis/streaklineDrawer/Particle.h"
#include "vis/streaklineDrawer/VelocityField.h"
#include "vis/streaklineDrawer/VelocitySiteData.h"

namespace hemelb
{
  namespace vis
  {
    namespace streaklinedrawer
    {
      class ParticleManager
      {
        public:
          ParticleManager(std::map<proc_t, NeighbouringProcessor>& iNeighbouringProcessors);

          void AddParticle(const Particle& iParticle);

          std::vector<Particle>& GetParticles();

          size_t GetNumberOfLocalParticles() const;

          void DeleteParticle(site_t iIndex);

          void DeleteAll();

          void ProcessParticleMovement();

          void CommunicateParticles(const geometry::LatticeData& iLatDat,
                                    proc_t iProcs,
                                    VelocityField& iVelocityField);

          /**
           * Function for debugging purposes. Logs information about the current state of the
           * particles.
           */
          void PrintParticleCount();

        private:
          std::vector<Particle> particles;
          std::map<proc_t, NeighbouringProcessor>& neighbouringProcessors;

      };
    }
  }
}

#endif // HEMELB_VIS_STREAKLINEDRAWER_PARTICLEMANAGER_H
