#include <cmath>

#include "debug/Debugger.h"
#include "vis/streaklineDrawer/ParticleManager.h"
#include "vis/streaklineDrawer/VelocityField.h"

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
          particles[i].x += particles[i].vx;
          particles[i].y += particles[i].vy;
          particles[i].z += particles[i].vz;
        }
      }

      // Communicate that particles current state to other processors.
      void ParticleManager::CommunicateParticles(const geometry::LatticeData& latticeData,
                                                 proc_t iProcs,
                                                 VelocityField& velocityField)
      {
        for (std::map<proc_t, NeighbouringProcessor>::iterator proc =
            neighbouringProcessors.begin(); proc != neighbouringProcessors.end(); proc++)
        {
          (*proc).second.PrepareToReceiveParticles();
        }

        unsigned int particles_temp = GetNumberOfLocalParticles();

        proc_t thisRank = topology::NetworkTopology::Instance()->GetLocalRank();

        for (int n = (int) (particles_temp - 1); n >= 0; n--)
        {
          site_t site_i = (site_t) particles[n].x;
          site_t site_j = (site_t) particles[n].y;
          site_t site_k = (site_t) particles[n].z;

          VelocitySiteData* siteVelocityData =
              velocityField.GetVelocitySiteData(latticeData, util::Vector3D<site_t>(site_i,
                                                                                    site_j,
                                                                                    site_k));

          if (siteVelocityData == NULL || siteVelocityData->proc_id == -1)
          {
            // TODO work out whether we should delete the particle here: DeleteParticle(n);
            continue;
          }
          else if (thisRank == siteVelocityData->proc_id)
          {
            continue;
          }

          neighbouringProcessors[siteVelocityData->proc_id].AddParticleToSend(particles[n]);
          DeleteParticle(n);
        }

        for (std::map<proc_t, NeighbouringProcessor>::iterator proc =
            neighbouringProcessors.begin(); proc != neighbouringProcessors.end(); proc++)
        {
          (*proc).second.PrepareToSendParticles();
        }

        for (std::map<proc_t, NeighbouringProcessor>::iterator proc =
            neighbouringProcessors.begin(); proc != neighbouringProcessors.end(); proc++)
        {
          (*proc).second.WaitForPreparationToReceiveParticles();
        }

        for (std::map<proc_t, NeighbouringProcessor>::iterator proc =
            neighbouringProcessors.begin(); proc != neighbouringProcessors.end(); proc++)
        {
          (*proc).second.SendParticles();
        }

        for (std::map<proc_t, NeighbouringProcessor>::iterator proc =
            neighbouringProcessors.begin(); proc != neighbouringProcessors.end(); proc++)
        {
          (*proc).second.WaitForParticlesToBeSent();
        }

        for (std::map<proc_t, NeighbouringProcessor>::iterator proc =
            neighbouringProcessors.begin(); proc != neighbouringProcessors.end(); proc++)
        {
          (*proc).second.ReceiveParticles();
        }

        for (std::map<proc_t, NeighbouringProcessor>::iterator proc =
            neighbouringProcessors.begin(); proc != neighbouringProcessors.end(); proc++)
        {
          NeighbouringProcessor& neighbourProc = (*proc).second;
          neighbourProc.WaitForParticlesToBeReceived();

          while (neighbourProc.ParticlesToBeRetrieved())
          {
            AddParticle(neighbourProc.PopNextReceivedParticle());
          }
        }
      }

      void ParticleManager::PrintParticleCount()
      {
        log::Logger::Log<log::Info, log::OnePerCore>("%i particles present", particles.size());
      }
    }
  }
}
