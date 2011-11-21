#ifndef HEMELB_VIS_STREAKLINEDRAWER_NEIGHBOURINGPROCESSOR_H
#define HEMELB_VIS_STREAKLINEDRAWER_NEIGHBOURINGPROCESSOR_H

#include <vector>

#include "mpiInclude.h"
#include "util/Vector3D.h"
#include "vis/streaklineDrawer/Particle.h"

namespace hemelb
{
  namespace vis
  {
    namespace streaklinedrawer
    {
      class NeighbouringProcessor
      {
        public:
          NeighbouringProcessor();

          NeighbouringProcessor(proc_t iID);

          void AddParticleToSend(const Particle& iParticle);

          bool ParticlesToBeRetrieved();

          const Particle& PopNextReceivedParticle();

          void PrepareToReceiveParticles();

          void PrepareToSendParticles();

          void WaitForPreparationToReceiveParticles();

          void SendParticles();

          void WaitForParticlesToBeSent();

          void ReceiveParticles();

          void WaitForParticlesToBeReceived();

          void AddSiteToRequestVelocityDataFor(site_t, site_t, site_t);

          void ExchangeSiteIds();

          void WaitForSiteIdExchange();

          void ExchangeVelocitiesForRequestedSites();

          void WaitForVelocityExchange();

          site_t GetSitesToReceiveCount();

          const util::Vector3D<float>& GetReceivedVelocity(const size_t receivedIndex) const;
          const util::Vector3D<site_t>& GetSiteCoordsBeingReceived(const size_t receivedIndex) const;

          void SetVelocityFieldToSend(const size_t sendIndex,
                                      const util::Vector3D<float>& velocityFieldToSend);

          size_t GetSendingSiteCount() const;

          util::Vector3D<site_t>& GetSiteCoorinates(size_t sendIndex);

          void ClearSendingSites();

        private:
          size_t numberOfParticlesToSend;
          std::vector<Particle> particlesToSend;

          size_t numberOfParticlesToReceive;
          std::vector<Particle> particlesToReceive;

          site_t receivedSitesCount;
          site_t sentSitesCount;

          std::vector<util::Vector3D<float> > sentVelocityField;
          std::vector<util::Vector3D<float> > receivedVelocityField;

          std::vector<util::Vector3D<site_t> > siteCoordsToSend;
          std::vector<util::Vector3D<site_t> > siteCoordsToReceive;

          proc_t neighbourRank;

          MPI_Request siteIdSendRequest;
          MPI_Request siteIdReceiveRequest;

          MPI_Request velocitySendRequest;
          MPI_Request velocityReceiveRequest;

          MPI_Request sendRequest;
          MPI_Request receiveRequest;
      };
    }
  }
}

#endif // HEMELB_VIS_STREAKLINEDRAWER_NEIGHBOURINGPROCESSOR_H
