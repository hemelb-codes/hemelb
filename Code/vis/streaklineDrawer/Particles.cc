#include <cmath>

#include "vis/streaklineDrawer/Particles.h"
#include "vis/streaklineDrawer/VelocityField.h"

namespace hemelb
{
  namespace vis
  {
    namespace streaklinedrawer
    {
      Particles::Particles(std::vector<NeighbouringProcessor>& iNeighbouringProcessors) :
          mNeighbouringProcessors(iNeighbouringProcessors)
      {
      }

      void Particles::AddParticle(const Particle& iParticle)
      {
        mParticles.push_back(iParticle);
      }

      std::vector<Particle>& Particles::GetParticles()
      {
        return mParticles;
      }

      size_t Particles::GetNumberOfParticles()
      {
        return mParticles.size();
      }

      void Particles::DeleteParticle(site_t iIndex)
      {
        assert(mParticles.size() > static_cast<size_t>(iIndex));

        //Move the particle at the end to position
        mParticles[iIndex] = mParticles.back();

        //Delete the now duplicated particle at the end
        mParticles.pop_back();
      }

      void Particles::DeleteAll()
      {
        mParticles.clear();
      }

      void Particles::ProcessParticleMovement()
      {
        for (unsigned int i = 0; i < GetNumberOfParticles(); i++)
        {
          // particle coords updating (dt = 1)
          mParticles[i].x += mParticles[i].vx;
          mParticles[i].y += mParticles[i].vy;
          mParticles[i].z += mParticles[i].vz;
        }
      }

      // Communicate that particles current state to other processors.
      void Particles::CommunicateParticles(const geometry::LatticeData& iLatDat,
                                           MPI_Request* iReq,
                                           proc_t iProcs,
                                           VelocityField& iVelocityField,
                                           proc_t* iFromProcIDToNeighbouringProcessorIndex)
      {
        for (size_t m = 0; m < mNeighbouringProcessors.size(); m++)
        {
	  mNeighbouringProcessors[m].PrepareToReceiveParticles();
	}
        /*for (size_t m = 0; m < mNeighbouringProcessors.size(); m++)
        {
          mNeighbouringProcessors[m].send_ps = 0;
	  }*/

        unsigned int particles_temp = GetNumberOfParticles();

        proc_t thisRank = topology::NetworkTopology::Instance()->GetLocalRank();

        for (int n = (int) (particles_temp - 1);n >= 0; n--) { 
	  site_t site_i = (unsigned int) mParticles[n].x;
          site_t site_j = (unsigned int) mParticles[n].y;
          site_t site_k = (unsigned int) mParticles[n].z;

          VelocitySiteData* vel_site_data_p = iVelocityField.velSiteDataPointer(iLatDat,
                                                                                site_i,
                                                                                site_j,
                                                                                site_k);

          if (vel_site_data_p == NULL || thisRank == vel_site_data_p->proc_id
              || vel_site_data_p->proc_id == -1)
          {
            continue;
          }
          proc_t m = iFromProcIDToNeighbouringProcessorIndex[vel_site_data_p->proc_id];

	  
          mNeighbouringProcessors[m].AddParticleToSend(mParticles[n]);
          DeleteParticle(n);
        }

        for (size_t m = 0; m < mNeighbouringProcessors.size(); m++)
        {
	  mNeighbouringProcessors[m].PrepareToSendParticles();
        }
	  
        for (size_t m = 0; m < mNeighbouringProcessors.size(); m++)
        {
          mNeighbouringProcessors[m].WaitForPreparationToReceiveParticles();
        }

        for (size_t m = 0; m < mNeighbouringProcessors.size(); m++)
        {
	  mNeighbouringProcessors[m].SendParticles();
        }
	for (size_t m = 0; m < mNeighbouringProcessors.size(); m++)
        {
	  mNeighbouringProcessors[m].WaitForParticlesToBeSent();
        }
        for (size_t m = 0; m < mNeighbouringProcessors.size(); m++)
        {
	  mNeighbouringProcessors[m].ReceiveParticles();
	}
	for (size_t m = 0; m < mNeighbouringProcessors.size(); m++)
        {
	  mNeighbouringProcessors[m].WaitForParticlesToBeReceived();

	  while (mNeighbouringProcessors[m].ParticlesToBeRetrieved())
	  {
	    AddParticle(mNeighbouringProcessors[m].RetrieveNextReceivedParticle());
	  }
	}
        
      }
    }
  }
}
