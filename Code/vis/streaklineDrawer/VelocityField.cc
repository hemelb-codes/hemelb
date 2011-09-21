#include <math.h>
#include <vector>
#include <cassert>

#include "vis/streaklineDrawer/StreaklineDrawer.h"
#include "vis/streaklineDrawer/VelocityField.h"
#include "vis/BlockTraverser.h"
#include "vis/SiteTraverser.h"

namespace hemelb
{
  namespace vis
  {
    namespace streaklinedrawer
    {
      VelocityField::VelocityField(std::vector<NeighProc>& iNeighProcs) :
	mNeighProcs(iNeighProcs)
      {
	counter = 1;
      }
      
      void VelocityField::BuildVelocityField(const geometry::LatticeData& iLatDat,
					     StreaklineDrawer* iStreaklineDrawer)
      {
	mVelocityField.resize(iLatDat.GetBlockCount());

	site_t inlet_sites = 0;

	BlockTraverser lBlockTraverser(iLatDat);
	do
	{
	  geometry::LatticeData::BlockData* lBlock = 
	    lBlockTraverser.GetCurrentBlockData();

	
	  if (lBlock->site_data == NULL)
	  {
	    continue;
	  }

	  SiteTraverser lSiteTraverser(iLatDat);
	  do
	  {
	    if (topology::NetworkTopology::Instance()->GetLocalRank() != 
		lBlock->ProcessorRankForEachBlockSite
		[lSiteTraverser.GetCurrentIndex()])
	    {
	      continue;
	    }

	    const site_t startI = util::NumericalFunctions::max<int>
	      (0, 
	       lBlockTraverser.GetX()*lBlockTraverser.GetBlockSize() +
	       lSiteTraverser.GetX() - 1);
	  
	    const site_t startJ = util::NumericalFunctions::max<int>
	      (0, 
	       lBlockTraverser.GetY()*lBlockTraverser.GetBlockSize() +
	       lSiteTraverser.GetY() - 1);
	
	    const site_t startK = util::NumericalFunctions::max<int>
	      (0, 
	       lBlockTraverser.GetZ()*lBlockTraverser.GetBlockSize() +
	       lSiteTraverser.GetZ() - 1);

	    const site_t endI =
	      util::NumericalFunctions::min<site_t>
	      (iLatDat.GetXSiteCount() - 1,
	       lBlockTraverser.GetX()*lBlockTraverser.GetBlockSize() +
	       lSiteTraverser.GetX() + 1);
	
	    const site_t endJ =
	      util::NumericalFunctions::min<site_t>
	      (iLatDat.GetYSiteCount() - 1,
	       lBlockTraverser.GetY()*lBlockTraverser.GetBlockSize() +
	       lSiteTraverser.GetY() + 1);
	
	    const site_t endK =
	      util::NumericalFunctions::min<site_t>
	      (iLatDat.GetZSiteCount() - 1,
	       lBlockTraverser.GetZ()*lBlockTraverser.GetBlockSize() +
	       lSiteTraverser.GetZ() + 1);
	

	    for (site_t neigh_i = startI; neigh_i <= endI; neigh_i++)
	    {
	      for (site_t neigh_j = startJ; neigh_j <= endJ; neigh_j++)
	      {
		for (site_t neigh_k = startK; neigh_k <= endK; neigh_k++)
		{
		  const proc_t* neigh_proc_id = iLatDat.GetProcIdFromGlobalCoords(neigh_i,
										  neigh_j,
										  neigh_k);

		  if (neigh_proc_id == NULL || *neigh_proc_id == BIG_NUMBER2)
		  {
		    continue;
		  }

		  initializeVelFieldBlock(iLatDat, neigh_i, neigh_j, neigh_k, *neigh_proc_id,
		    iStreaklineDrawer);

		  if (topology::NetworkTopology::Instance()->GetLocalRank() == *neigh_proc_id)
		  {
		    continue;
		  }

		  VelocitySiteData* vel_site_data_p = velSiteDataPointer(iLatDat,
									 neigh_i,
									 neigh_j,
									 neigh_k);

		  if (vel_site_data_p->counter == counter)
		  {
		    continue;
		  }

		  vel_site_data_p->counter = counter;

		  bool seenSelf = false;
		  for (size_t mm = 0; mm < mNeighProcs.size() && !seenSelf; mm++)
		  {
		    if (*neigh_proc_id == mNeighProcs[mm].id)
		    {
		      seenSelf = true;
		      ++mNeighProcs[mm].send_vs;
		    }
		  }
		  if (seenSelf)
		  {
		    continue;
		  }

		  NeighProc lNew;

		  lNew.id = *neigh_proc_id;
		  lNew.send_vs = 1;
		  mNeighProcs.push_back(lNew);
		}
	      }
	    }

	    site_t lSiteIndex = lBlock->site_data[lSiteTraverser.GetCurrentIndex()];

	    // if the lattice site is not an inlet
	    if (iLatDat.GetSiteType(lSiteIndex) != geometry::LatticeData::INLET_TYPE)
	    {
	      continue;
	    }
	    ++inlet_sites;

	    if (inlet_sites % 50 != 0)
	    {
	      continue;
	    }

	    iStreaklineDrawer->mParticleSeeds.push_back(
	      Particle( static_cast<float>(lBlockTraverser.GetX()*lBlockTraverser.GetBlockSize() +
					   lSiteTraverser.GetX()),
			static_cast<float>(lBlockTraverser.GetY()*lBlockTraverser.GetBlockSize() +
					   lSiteTraverser.GetY()),
			static_cast<float>(lBlockTraverser.GetZ()*lBlockTraverser.GetBlockSize() +
					   lSiteTraverser.GetZ()), 
			iLatDat.GetBoundaryId(lSiteIndex) )
	      );
	    

	  }
	  while (lSiteTraverser.TraverseOne());
	}
	while (lBlockTraverser.TraverseOne());

	for (site_t n = 0; n < iLatDat.GetBlockCount(); n++)
	{
	  if (mVelocityField[n].empty())
	  {
	    continue;
	  }
	
	  for (site_t m = 0; m < iLatDat.GetSitesPerBlockVolumeUnit(); m++)
	  {
	    mVelocityField[n][m].counter = counter;
	  }
	  if (iLatDat.GetBlock(n)->site_data == NULL)
	    continue;

	  for (site_t m = 0; m < iLatDat.GetSitesPerBlockVolumeUnit(); m++)
	  {
	    mVelocityField[n][m].site_id = iLatDat.GetBlock(n)->site_data[m];
	  }
	}

	counter = 0;
      }
     
      
      bool VelocityField::BlockContainsData(size_t iBlockNumber)
      {
	return !mVelocityField[iBlockNumber].empty();
      }

      VelocitySiteData& VelocityField::GetSiteData(size_t iBlockNumber, size_t iSiteNumber)
      {
	return mVelocityField[iBlockNumber][iSiteNumber];
      }
      
      // Returns the velocity site data for a given index, or NULL if the index isn't valid / has
      // no data.
       VelocitySiteData* VelocityField::velSiteDataPointer
       ( const geometry::LatticeData& iLatDat,
	 site_t site_i,
	 site_t site_j,
	 site_t site_k )
       {
	if (site_i >= iLatDat.GetXSiteCount() || site_j >= iLatDat.GetYSiteCount() || site_k
	    >= iLatDat.GetZSiteCount())
	{
	  return NULL;
	}
	site_t i = site_i >> iLatDat.GetLog2BlockSize();
	site_t j = site_j >> iLatDat.GetLog2BlockSize();
	site_t k = site_k >> iLatDat.GetLog2BlockSize();

	site_t block_id = iLatDat.GetBlockIdFromBlockCoords(i, j, k);

	if (!BlockContainsData(static_cast<size_t>(block_id)))
	{
	  return NULL;
	}
	site_t ii = site_i - (i << iLatDat.GetLog2BlockSize());
	site_t jj = site_j - (j << iLatDat.GetLog2BlockSize());
	site_t kk = site_k - (k << iLatDat.GetLog2BlockSize());

	site_t site_id =
	  ( ( (ii << iLatDat.GetLog2BlockSize()) + jj) << iLatDat.GetLog2BlockSize()) + kk;

	return &GetSiteData(block_id,site_id);
      }

      // Function to initialise the velocity field at given coordinates.
      void VelocityField::initializeVelFieldBlock(const geometry::LatticeData& iLatDat,
						  site_t site_i,
						  site_t site_j,
						  site_t site_k,
						  proc_t proc_id,
						  StreaklineDrawer* iStreaklineDrawer)
      {
	site_t i = site_i >> iLatDat.GetLog2BlockSize();
	site_t j = site_j >> iLatDat.GetLog2BlockSize();
	site_t k = site_k >> iLatDat.GetLog2BlockSize();

	site_t block_id = iLatDat.GetBlockIdFromBlockCoords(i, j, k);

	if (!BlockContainsData(block_id))
	{
	  mVelocityField[block_id] = 
	    std::vector<VelocitySiteData>(iLatDat.GetSitesPerBlockVolumeUnit());
	}

	site_t ii = site_i - (i << iLatDat.GetLog2BlockSize());
	site_t jj = site_j - (j << iLatDat.GetLog2BlockSize());
	site_t kk = site_k - (k << iLatDat.GetLog2BlockSize());

	site_t site_id =
	  ( ( (ii << iLatDat.GetLog2BlockSize()) + jj) << iLatDat.GetLog2BlockSize()) + kk;
	mVelocityField[block_id][site_id].proc_id = proc_id;
      }
    
      // Populate the matrix v with all the velocity field data at each index.
      void VelocityField::localVelField(site_t iX,
					site_t iY,
					site_t iZ,
					float v[2][2][2][3],
					int *is_interior,
					const geometry::LatticeData& iLatDat,
					StreaklineDrawer* iStreaklineDrawer)
      {
	/*
	site_t site_i = (site_t) iStreaklineDrawer->mParticleVec[p_index].x;
	site_t site_j = (site_t) iStreaklineDrawer->mParticleVec[p_index].y;
	site_t site_k = (site_t) iStreaklineDrawer->mParticleVec[p_index].z;*/

	*is_interior = 1;

	const proc_t thisRank = topology::NetworkTopology::Instance()->GetLocalRank();

	for (unsigned int i = 0; i < 2; i++)
	{
	  site_t neigh_i = iX + i;

	  for (unsigned int j = 0; j < 2; j++)
	  {
	    site_t neigh_j = iY + j;

	    for (unsigned int k = 0; k < 2; k++)
	    {
	      site_t neigh_k = iZ + k;

	      VelocitySiteData *vel_site_data_p = 
		velSiteDataPointer(iLatDat, neigh_i, neigh_j, neigh_k);

	      if (vel_site_data_p == NULL || vel_site_data_p->proc_id == -1)
	      {
		// it is a solid site and the velocity is
		// assumed to be zero
		v[i][j][k][0] = v[i][j][k][1] = v[i][j][k][2] = 0.0F;
		continue;
	      }
	      if (thisRank != vel_site_data_p->proc_id)
	      {
		*is_interior = 0;
	      }
	      if (vel_site_data_p->counter == counter)
	      {
		// This means that the local velocity has already been
		// calculated at the current time step if the site
		// belongs to the current processor; if not, the
		// following instructions have no effect
		v[i][j][k][0] = vel_site_data_p->vx;
		v[i][j][k][1] = vel_site_data_p->vy;
		v[i][j][k][2] = vel_site_data_p->vz;
	      }
	      else if (thisRank == vel_site_data_p->proc_id)
	      {
		// the local counter is set equal to the global one
		// and the local velocity is calculated
		vel_site_data_p->counter = counter;
		distribn_t density, vx, vy, vz;

		D3Q15::CalculateDensityAndVelocity(iLatDat.GetFOld(vel_site_data_p->site_id
								   * D3Q15::NUMVECTORS), density, vx, vy, vz);

		v[i][j][k][0] = vel_site_data_p->vx = (float) (vx / density);
		v[i][j][k][1] = vel_site_data_p->vy = (float) (vy / density);
		v[i][j][k][2] = vel_site_data_p->vz = (float) (vz / density);
	      }
	      else
	      {
		vel_site_data_p->counter = counter;

		proc_t m = iStreaklineDrawer->from_proc_id_to_neigh_proc_index[vel_site_data_p->proc_id];

		mNeighProcs[m].s_to_send[3 * mNeighProcs[m].send_vs + 0] = neigh_i;
		mNeighProcs[m].s_to_send[3 * mNeighProcs[m].send_vs + 1] = neigh_j;
		mNeighProcs[m].s_to_send[3 * mNeighProcs[m].send_vs + 2] = neigh_k;
		++(mNeighProcs[m].send_vs);
	      }
	    }
	  }
	}
      }

     // Interpolates a velocity field to get the velocity at the position of a particle.
      void VelocityField::GetVelocityAtPoint(float x,
					     float y,
					     float z,
					     float v[2][2][2][3],
					     float interp_v[3])
      {
	float v_00z, v_01z, v_10z, v_11z, v_0y, v_1y;

	float dummy;

	float dx = modff(x, &dummy);
	float dy = modff(y, &dummy);
	float dz = modff(z, &dummy);

	for (int l = 0; l < 3; l++)
	{
	  v_00z = (1.F - dz) * v[0][0][0][l] + dz * v[0][0][1][l];
	  v_01z = (1.F - dz) * v[0][1][0][l] + dz * v[0][1][1][l];
	  v_10z = (1.F - dz) * v[1][0][0][l] + dz * v[1][0][1][l];
	  v_11z = (1.F - dz) * v[1][1][0][l] + dz * v[1][1][1][l];

	  v_0y = (1.F - dy) * v_00z + dy * v_01z;
	  v_1y = (1.F - dy) * v_10z + dy * v_11z;

	  interp_v[l] = (1.F - dx) * v_0y + dx * v_1y;
	}
      }

    }
  }
}
