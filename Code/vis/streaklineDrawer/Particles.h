#ifndef HEMELB_VIS_PARTICLES_H
#define HEMELB_VIS_PARTICLES_H

#include <vector>

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
      class Particles
      {
        public:
          Particles(std::vector<NeighbouringProcessor>& iNeighbouringProcessors);

          void AddParticle(const Particle& iParticle);

          std::vector<Particle>& GetParticles();

          size_t GetNumberOfParticles();

          void DeleteParticle(site_t iIndex);

          void DeleteAll();

          void ProcessParticleMovement();

          void CommunicateParticles(const geometry::LatticeData& iLatDat,
                                    MPI_Request* iReq,
                                    proc_t iProcs,
                                    VelocityField& iVelocityField,
                                    proc_t* iFromProcIDToNeighbouringProcessorIndex);

        private:
          std::vector<Particle> mParticles;
          std::vector<NeighbouringProcessor>& mNeighbouringProcessors;

      };
    }
  }
}

#endif // HEMELB_VIS_STREAKLINEDRAWER_H
