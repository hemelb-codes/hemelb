// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include <cmath>

#include "vis/streaklineDrawer/ParticleManager.h"

namespace hemelb
{
  namespace vis
  {
    namespace streaklinedrawer
    {
      ParticleManager::ParticleManager(std::map<proc_t, NeighbouringProcessor>& iNeighbouringProcessors) :
        neighbouringProcessors(iNeighbouringProcessors)
      {
      }

      void ParticleManager::AddParticle(const Particle& iParticle)
      {
        particles.push_back(iParticle);
      }

      std::vector<Particle>& ParticleManager::GetParticles()
      {
        return particles;
      }

      size_t ParticleManager::GetNumberOfLocalParticles() const
      {
        return particles.size();
      }

      void ParticleManager::DeleteParticle(site_t iIndex)
      {
        assert(particles.size() > static_cast<size_t> (iIndex));

        //Move the particle at the end to position
        particles[iIndex] = particles.back();

        //Delete the now duplicated particle at the end
        particles.pop_back();
      }

      void ParticleManager::DeleteAll()
      {
        particles.clear();
      }

      void ParticleManager::ProcessParticleMovement()
      {
        for (unsigned int i = 0; i < GetNumberOfLocalParticles(); i++)
        {
          // particle coords updating (dt = 1)
          particles[i].position += particles[i].velocity;
        }
      }

      // Communicate the particles' current state to other processors.
      void ParticleManager::CommunicateParticles(net::Net& streakNet,
                                                 const geometry::LatticeData& latticeData,
                                                 VelocityField& velocityField)
      {
        unsigned int particles_temp = GetNumberOfLocalParticles();

        proc_t thisRank = net::IOCommunicator::Instance()->Rank();

        for (int n = (int) (particles_temp - 1); n >= 0; n--)
        {
          VelocitySiteData* siteVelocityData =
              velocityField.GetVelocitySiteData(latticeData,
                                                util::Vector3D<site_t>(particles[n].position));

          // TODO can we get rid of the first test?
          if (siteVelocityData == NULL || siteVelocityData->proc_id == -1)
          {
            continue;
          }
          else if (thisRank == siteVelocityData->proc_id)
          {
            continue;
          }

          neighbouringProcessors[siteVelocityData->proc_id].AddParticleToSend(particles[n]);
          DeleteParticle(n);
        }

        net::Net net;

        for (std::map<proc_t, NeighbouringProcessor>::iterator proc =
            neighbouringProcessors.begin(); proc != neighbouringProcessors.end(); ++proc)
        {
          (*proc).second.ExchangeParticleCounts(streakNet);
        }
        streakNet.Dispatch();

        for (std::map<proc_t, NeighbouringProcessor>::iterator proc =
            neighbouringProcessors.begin(); proc != neighbouringProcessors.end(); ++proc)
        {
          (*proc).second.ExchangeParticles(streakNet);
        }

        streakNet.Dispatch();

        for (std::map<proc_t, NeighbouringProcessor>::iterator proc =
            neighbouringProcessors.begin(); proc != neighbouringProcessors.end(); ++proc)
        {
          NeighbouringProcessor& neighbourProc = (*proc).second;
          neighbourProc.ClearParticleSendingList();

          while (neighbourProc.ParticlesToBeRetrieved())
          {
            AddParticle(neighbourProc.PopNextReceivedParticle());
          }
        }
      }

    }
  }
}
