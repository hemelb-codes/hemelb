// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "geometry/neighbouring/NeighbouringDomain.h"
#include "geometry/neighbouring/NeighbouringSite.h"
#include "log/Logger.h"
namespace hemelb
{
  namespace geometry
  {
    namespace neighbouring
    {

      NeighbouringDomain::NeighbouringDomain(const lb::lattices::LatticeInfo& latticeInfo) :
          distanceToWall(), wallNormalAtSite(), siteData(),
              latticeInfo(latticeInfo)
      {
      }

      void NeighbouringDomain::SaveSite(site_t index,
                                        const std::vector<distribn_t> &distances,
                                        const util::Vector3D<distribn_t> &normal,
                                        const SiteData & data)
      {
        for (unsigned int direction = 0; direction < latticeInfo.GetNumVectors() - 1; direction++)
        {
          GetCutDistances(index)[direction] = distances[direction];
        }
        GetNormalToWall(index) = normal;
        GetSiteData(index) = data;
      }

      auto NeighbouringDomain::GetSite(site_t globalIndex) -> NeighbouringSite
      {
        return NeighbouringSite(globalIndex, *this);
      }

      auto NeighbouringDomain::GetSite(site_t globalIndex) const -> ConstNeighbouringSite
      {
        return ConstNeighbouringSite(globalIndex, *this);
      }

      const util::Vector3D<distribn_t>& NeighbouringDomain::GetNormalToWall(
          site_t globalIndex) const
      {
        return wallNormalAtSite.find(globalIndex)->second;
      }

      util::Vector3D<distribn_t>& NeighbouringDomain::GetNormalToWall(site_t globalIndex)
      {
        return wallNormalAtSite[globalIndex];
      }

      auto NeighbouringFieldData::GetSite(site_t globalIndex) -> NeighbouringSite
      {
          return NeighbouringSite(globalIndex, *this);
      }

      auto NeighbouringFieldData::GetSite(site_t globalIndex) const -> ConstNeighbouringSite
      {
          return ConstNeighbouringSite(globalIndex, *this);
      }

      distribn_t* NeighbouringFieldData::GetFOld(site_t distributionIndex)
      {
        const auto Q = m_domain->latticeInfo.GetNumVectors();
        site_t globalIndex = distributionIndex / Q;
        site_t direction = distributionIndex % Q;
        auto& buffer = distributions[globalIndex];
        buffer.resize(Q);
        return &buffer[direction];
      }

      std::vector<distribn_t>& NeighbouringFieldData::GetDistribution(site_t globalIndex)
      {
        return distributions[globalIndex];
      }

      const distribn_t* NeighbouringFieldData::GetFOld(site_t distributionIndex) const
      {
        const auto Q = m_domain->latticeInfo.GetNumVectors();
        site_t globalIndex = distributionIndex / Q;
        site_t direction = distributionIndex % Q;
        return &distributions.find(globalIndex)->second[direction];
      }

      const SiteData & NeighbouringDomain::GetSiteData(site_t globalIndex) const
      {
        return siteData.find(globalIndex)->second;
      }

      SiteData& NeighbouringDomain::GetSiteData(site_t globalIndex)
      {
        return siteData[globalIndex];
      }

      const distribn_t * NeighbouringDomain::GetCutDistances(site_t globalIndex) const
      {
        return &distanceToWall.find(globalIndex)->second.front();
      }

      distribn_t* NeighbouringDomain::GetCutDistances(site_t globalIndex)
      {
        std::vector<distribn_t> &buffer = distanceToWall[globalIndex];
        buffer.resize(latticeInfo.GetNumVectors() - 1);
        return &buffer.front();
      }

        void NeighbouringFieldData::SaveSite(site_t index,
                                          const std::vector<distribn_t> &distribution,
                                          const std::vector<distribn_t> &distances,
                                          const util::Vector3D<distribn_t> &normal,
                                          const SiteData & data) {
            GetDistribution(index) = distribution;
            m_domain->SaveSite(index, distances, normal, data);
        }
    }
  }
}

