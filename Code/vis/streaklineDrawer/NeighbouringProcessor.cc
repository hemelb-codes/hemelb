// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include <vector>
#include <cassert>

#include "constants.h"
#include "vis/streaklineDrawer/NeighbouringProcessor.h"

namespace hemelb
{
  namespace vis
  {
    namespace streaklinedrawer
    {

      NeighbouringProcessor::NeighbouringProcessor()
      {
      }

      NeighbouringProcessor::NeighbouringProcessor(proc_t neighbourRankIn) :
          neighbourRank(neighbourRankIn)
      {
      }

      void NeighbouringProcessor::AddParticleToSend(const Particle& particle)
      {
        particlesToSend.push_back(particle);
      }

      bool NeighbouringProcessor::ParticlesToBeRetrieved()
      {
        return (particlesToReceive.size() > 0);
      }

      const Particle& NeighbouringProcessor::PopNextReceivedParticle()
      {
        const Particle& lParticle = particlesToReceive.back();
        particlesToReceive.pop_back();
        return lParticle;
      }

      void NeighbouringProcessor::ExchangeParticleCounts(net::Net& net)
      {
        numberOfParticlesToSend = particlesToSend.size();
        net.RequestSendR(numberOfParticlesToSend, neighbourRank);
        net.RequestReceiveR(numberOfParticlesToReceive, neighbourRank);
      }

      void NeighbouringProcessor::ClearParticleSendingList()
      {
        particlesToSend.clear();
      }

      void NeighbouringProcessor::ExchangeParticles(net::Net& net)
      {
        if (numberOfParticlesToReceive > 0)
        {
          particlesToReceive.resize(numberOfParticlesToReceive);

          net.RequestReceiveV(particlesToReceive, neighbourRank);
        }

        if (particlesToSend.size() > 0)
        {
          net.RequestSendV(particlesToSend, neighbourRank); //Request
        }
      }

      void NeighbouringProcessor::AddSiteToRequestVelocityDataFor(site_t siteI, site_t siteJ,
                                                                  site_t siteK)
      {
        siteCoordsRequestedByThisCore.push_back(util::Vector3D<site_t>(siteI, siteJ, siteK));
      }

      void NeighbouringProcessor::ExchangeSiteIdCounts(net::Net& net)
      {
        numberOfSitesRequestedByThisCore = siteCoordsRequestedByThisCore.size();

        net.RequestReceiveR(numberOfSiteBeingRequestedByNeighbour, neighbourRank);
        net.RequestSendR(numberOfSitesRequestedByThisCore, neighbourRank);
      }

      void NeighbouringProcessor::ExchangeSiteIds(net::Net& net)
      {
        if (numberOfSiteBeingRequestedByNeighbour > 0)
        {
          siteCoordsRequestedByNeighbour.resize(numberOfSiteBeingRequestedByNeighbour);
          velocityFieldDataForNeighbour.resize(numberOfSiteBeingRequestedByNeighbour);

          net.RequestReceive(&siteCoordsRequestedByNeighbour[0],
                             (int) numberOfSiteBeingRequestedByNeighbour,
                             neighbourRank);
        }

        if (numberOfSitesRequestedByThisCore > 0)
        {
          velocityFieldDataFromNeighbour.resize(numberOfSitesRequestedByThisCore);

          net.RequestSendV(siteCoordsRequestedByThisCore, neighbourRank);
        }
      }

      void NeighbouringProcessor::ExchangeVelocitiesForRequestedSites(net::Net& net)
      {
        if (numberOfSitesRequestedByThisCore > 0)
        {
          velocityFieldDataFromNeighbour.resize(numberOfSitesRequestedByThisCore);

          net.RequestReceiveV(velocityFieldDataFromNeighbour, neighbourRank);
        }

        if (numberOfSiteBeingRequestedByNeighbour > 0)
        {
          velocityFieldDataForNeighbour.resize(numberOfSiteBeingRequestedByNeighbour);

          net.RequestSendV(velocityFieldDataForNeighbour, neighbourRank);
        }
      }

      site_t NeighbouringProcessor::GetNumberOfSitesRequestedByNeighbour() const
      {
        return siteCoordsRequestedByNeighbour.size();
      }

      const util::Vector3D<float>& NeighbouringProcessor::GetReceivedVelocityField(
          const site_t receivedIndex) const
      {
        return velocityFieldDataFromNeighbour[receivedIndex];
      }

      const util::Vector3D<site_t>& NeighbouringProcessor::GetSiteCoordsBeingRequestedByNeighbour(
          const site_t receivedIndex) const
      {
        return siteCoordsRequestedByNeighbour[receivedIndex];
      }

      void NeighbouringProcessor::SetVelocityFieldToSend(
          const site_t sendIndex, const util::Vector3D<float>& velocityFieldToSend)
      {
        velocityFieldDataForNeighbour[sendIndex] = velocityFieldToSend;
      }

      site_t NeighbouringProcessor::GetNumberOfSitesRequestedByThisCore() const
      {
        return siteCoordsRequestedByThisCore.size();
      }

      const util::Vector3D<site_t>& NeighbouringProcessor::GetSendingSiteCoorinates(
          site_t sendIndex) const
      {
        return siteCoordsRequestedByThisCore[sendIndex];
      }

      void NeighbouringProcessor::ClearListOfRequestedSites()
      {
        siteCoordsRequestedByThisCore.clear();
      }

    }
  }
}
