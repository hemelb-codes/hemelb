#include "NeighbouringLatticeData.h"
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
        SaveDistribution(index, distribution);
        SaveDistances(index, distances);
        SaveNormal(index, normal);
        SaveData(index, data);
      }

      void NeighbouringLatticeData::SaveDistribution(site_t index, const std::vector<distribn_t> &distribution)
      {
        distributions[index] = distribution;
      }

      void NeighbouringLatticeData::SaveDistances(site_t index, const std::vector<distribn_t> & distances)
      {
        distanceToWall[index] = distances;
      }

      void NeighbouringLatticeData::SaveNormal(site_t index, const util::Vector3D<distribn_t> &normal)
      {
        wallNormalAtSite[index] = normal;
      }

      void NeighbouringLatticeData::SaveData(site_t index, const SiteData & data)
      {
        // no default constructor for site data, so don't use operator[]
        siteData.insert(std::pair<site_t, SiteData>(index, data));
      }


      NeighbouringSite NeighbouringLatticeData::GetSite(site_t globalIndex)
      {
        return NeighbouringSite(globalIndex, *this);
      }

      const util::Vector3D<distribn_t>& NeighbouringLatticeData::GetNormalToWall(site_t globalIndex) const
      {
        return wallNormalAtSite.find(globalIndex)->second;
      }

      distribn_t* NeighbouringLatticeData::GetFOld(site_t distributionIndex)
      {
        site_t globalIndex = distributionIndex / latticeInfo.GetNumVectors();
        site_t direction = distributionIndex % latticeInfo.GetNumVectors();
        return &distributions[globalIndex][direction];
      }
      const distribn_t* NeighbouringLatticeData::GetFOld(site_t distributionIndex) const
      {
        site_t globalIndex = distributionIndex / latticeInfo.GetNumVectors();
        site_t direction = distributionIndex % latticeInfo.GetNumVectors();
        return &distributions.find(globalIndex)->second[direction];
      }

      SiteData NeighbouringLatticeData::GetSiteData(site_t globalIndex) const
      {
        return siteData.find(globalIndex)->second;
      }

    }
  }
}

