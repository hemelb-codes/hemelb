#include <cmath>

#include "vis/streaklineDrawer/Particles.h"
#include "vis/streaklineDrawer/VelocityField.h"

namespace hemelb
{
  namespace vis
  {
    namespace streaklinedrawer
    {
      Particles::Particles(std::vector<NeighProc>& iNeighProcs) :
      	mNeighProcs(iNeighProcs)
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
	assert(mParticles.size() > iIndex);

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
	for (unsigned int i= 0; i < GetNumberOfParticles(); i++)
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
					   proc_t* iFromProcIDToNeighProcIndex)
      {
	for (size_t m = 0; m < mNeighProcs.size(); m++)
	{
	  MPI_Irecv(&mNeighProcs[m].recv_ps,
		    1,
		    MpiDataType<site_t> (),
		    mNeighProcs[m].id,
		    30,
		    MPI_COMM_WORLD,
		    &iReq[iProcs + mNeighProcs[m].id]);
	}
	for (size_t m = 0; m < mNeighProcs.size(); m++)
	{
	  mNeighProcs[m].send_ps = 0;
	}

	unsigned int particles_temp = GetNumberOfParticles();

	proc_t thisRank = topology::NetworkTopology::Instance()->GetLocalRank();

	for (int n = (int) (particles_temp - 1); n >= 0; n--)
	{
	  site_t site_i = (unsigned int) mParticles[n].x;
	  site_t site_j = (unsigned int) mParticles[n].y;
	  site_t site_k = (unsigned int) mParticles[n].z;

	  VelocitySiteData* vel_site_data_p = 
	    iVelocityField.velSiteDataPointer(iLatDat, site_i, site_j, site_k);

	  if (vel_site_data_p == NULL || thisRank == vel_site_data_p->proc_id
	      || vel_site_data_p->proc_id == -1)
	  {
	    continue;
	  }
	  proc_t m = iFromProcIDToNeighProcIndex[vel_site_data_p->proc_id];


	  mNeighProcs[m].p_to_send[5 * mNeighProcs[m].send_ps + 0] = mParticles[n].x;
	  mNeighProcs[m].p_to_send[5 * mNeighProcs[m].send_ps + 1] = mParticles[n].y;
	  mNeighProcs[m].p_to_send[5 * mNeighProcs[m].send_ps + 2] = mParticles[n].z;
	  mNeighProcs[m].p_to_send[5 * mNeighProcs[m].send_ps + 3] = mParticles[n].vel;
	  mNeighProcs[m].p_to_send[5 * mNeighProcs[m].send_ps + 4] = (float) mParticles[n].inletID
	    + 0.1F;
	  ++mNeighProcs[m].send_ps;

	  DeleteParticle(n);
	}
	for (size_t m = 0; m < mNeighProcs.size(); m++)
	{
	  MPI_Isend(&mNeighProcs[m].send_ps,
		    1,
		    MpiDataType<site_t> (),
		    mNeighProcs[m].id,
		    30,
		    MPI_COMM_WORLD,
		    &iReq[mNeighProcs[m].id]);

	}
	for (size_t m = 0; m < mNeighProcs.size(); m++)
	{
	  MPI_Wait(&iReq[iProcs + mNeighProcs[m].id], MPI_STATUS_IGNORE);
	}

	for (size_t m = 0; m < mNeighProcs.size(); m++)
	{
	  if (mNeighProcs[m].send_ps > 0)
	  {
	    MPI_Isend(&mNeighProcs[m].p_to_send[0],
		      (int) mNeighProcs[m].send_ps * 5,
		      MpiDataType(mNeighProcs[m].p_to_send[0]),
		      mNeighProcs[m].id,
		      40,
		      MPI_COMM_WORLD,
		      &iReq[mNeighProcs[m].id]);

	    MPI_Wait(&iReq[mNeighProcs[m].id], MPI_STATUS_IGNORE);
	  }
	}
	for (size_t m = 0; m < mNeighProcs.size(); m++)
	{
	  if (mNeighProcs[m].recv_ps > 0)
	  {
	    MPI_Irecv(&mNeighProcs[m].p_to_recv[0],
		      (int) mNeighProcs[m].recv_ps * 5,
		      MpiDataType(mNeighProcs[m].p_to_recv[0]),
		      mNeighProcs[m].id,
		      40,
		      MPI_COMM_WORLD,
		      &iReq[iProcs + mNeighProcs[m].id]);
	    MPI_Wait(&iReq[iProcs + mNeighProcs[m].id], MPI_STATUS_IGNORE);

	    for (proc_t n = 0; n < mNeighProcs[m].recv_ps; n++)
	    {
	      AddParticle(Particle(mNeighProcs[m].p_to_recv[5 * n + 0],
				   mNeighProcs[m].p_to_recv[5 * n + 1],
				   mNeighProcs[m].p_to_recv[5 * n + 2],
				   mNeighProcs[m].p_to_recv[5 * n + 3],
				   (int) mNeighProcs[m].p_to_recv[5 * n + 4]));
	    }
	  }
	}
      }
    }
  }
}
