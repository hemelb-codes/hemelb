
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "geometry/neighbouring/NeighbouringLatticeData.h"
#include "geometry/neighbouring/NeighbouringSite.h"
#include "log/Logger.h"
namespace hemelb
{
  namespace geometry
  {
    namespace neighbouring
    {

      NeighbouringLatticeData::NeighbouringLatticeData(const lb::lattices::LatticeInfo& latticeInfo) :
          distributions(), distanceToWall(), wallNormalAtSite(), siteData(), latticeInfo(latticeInfo)
      {
      }

      void NeighbouringLatticeData::SaveSite(site_t index,
                                             const std::vector<distribn_t> &distribution,
                                             const std::vector<distribn_t> &distances,
                                             const util::Vector3D<distribn_t> &normal,
                                             const SiteData & data)
      {
        GetDistribution(index) = distribution;
        for (unsigned int direction = 0; direction < latticeInfo.GetNumVectors() - 1; direction++)
        {
          GetCutDistances(index)[direction] = distances[direction];
        }
        GetNormalToWall(index) = normal;
        GetSiteData(index) = data;
      }

      NeighbouringSite NeighbouringLatticeData::GetSite(site_t globalIndex)
      {
        return NeighbouringSite(globalIndex, *this);
      }

      ConstNeighbouringSite NeighbouringLatticeData::GetSite(site_t globalIndex) const
      {
        return ConstNeighbouringSite(globalIndex, *this);
      }

      const util::Vector3D<distribn_t>& NeighbouringLatticeData::GetNormalToWall(site_t globalIndex) const
      {
        return wallNormalAtSite.find(globalIndex)->second;
      }

      util::Vector3D<distribn_t>& NeighbouringLatticeData::GetNormalToWall(site_t globalIndex)
      {
        return wallNormalAtSite[globalIndex];
      }

      distribn_t* NeighbouringLatticeData::GetFOld(site_t distributionIndex)
      {
        site_t globalIndex = distributionIndex / latticeInfo.GetNumVectors();
        site_t direction = distributionIndex % latticeInfo.GetNumVectors();
        std::vector<distribn_t> &buffer = distributions[globalIndex];
        buffer.resize(latticeInfo.GetNumVectors());
        return &buffer[direction];
      }

      std::vector<distribn_t>& NeighbouringLatticeData::GetDistribution(site_t globalIndex)
      {
        return distributions[globalIndex];
      }

      const distribn_t* NeighbouringLatticeData::GetFOld(site_t distributionIndex) const
      {
        site_t globalIndex = distributionIndex / latticeInfo.GetNumVectors();
        site_t direction = distributionIndex % latticeInfo.GetNumVectors();
        return &distributions.find(globalIndex)->second[direction];
      }

      const SiteData & NeighbouringLatticeData::GetSiteData(site_t globalIndex) const
      {
        return siteData.find(globalIndex)->second;
      }

      SiteData& NeighbouringLatticeData::GetSiteData(site_t globalIndex)
      {
        return siteData[globalIndex];
      }

      const distribn_t * NeighbouringLatticeData::GetCutDistances(site_t globalIndex) const
      {
        return &distanceToWall.find(globalIndex)->second.front();
      }

      distribn_t* NeighbouringLatticeData::GetCutDistances(site_t globalIndex)
      {
        std::vector<distribn_t> &buffer = distanceToWall[globalIndex];
        buffer.resize(latticeInfo.GetNumVectors() - 1);
        return &buffer.front();
      }

    }
  }
}

