#ifndef HEMELB_GEOMETRY_NEIGHBOURING_NEIGHBOURINGLATTICEDATA_H
#define HEMELB_GEOMETRY_NEIGHBOURING_NEIGHBOURINGLATTICEDATA_H
#include <map>
#include "geometry/SiteData.h"
#include "geometry/neighbouring/NeighbouringSite.h"
namespace hemelb
{
  namespace geometry
  {
    namespace neighbouring
    {
      // Here, all site indices are GLOBAL index.
      // Local users must determine the global index of the site they are interested in.

      class NeighbouringLatticeData
      {

        public:
          friend class BaseSite<NeighbouringLatticeData> ; //! Let the inner classes have access to site-related data that's otherwise private.

          NeighbouringLatticeData(const lb::lattices::LatticeInfo& latticeInfo) :
              distributions(), distanceToWall(), wallNormalAtSite(), siteData(), latticeInfo(latticeInfo)
          {
          }

          virtual ~NeighbouringLatticeData()
          {
          }

          void SaveSite(site_t index,
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
          void SaveDistribution(site_t index, const std::vector<distribn_t> &distribution)
          {
            distributions[index] = distribution;
          }

          void SaveDistances(site_t index, const std::vector<distribn_t> & distances)
          {
            distanceToWall[index] = distances;
          }

          void SaveNormal(site_t index, const util::Vector3D<distribn_t> &normal)
          {
            wallNormalAtSite[index] = normal;
          }

          void SaveData(site_t index, const SiteData & data)
          {
            // no default constructor for site data, so don't use operator[]
            siteData.insert(std::pair<site_t, SiteData>(index, data));
          }

          /**
           * Get a site object for the given index.
           * @param localIndex
           * @return
           */
          inline NeighbouringSite GetSite(site_t globalIndex)
          {
            return NeighbouringSite(globalIndex, *this);
          }

          /**
           * Get the wall normal at the given site
           * @param iSiteIndex
           * @return
           */
          inline const util::Vector3D<distribn_t>& GetNormalToWall(site_t globalIndex) const
          {
            return wallNormalAtSite.find(globalIndex)->second;
          }

          /**
           * Get a pointer to the fOld array starting at the requested index
           * LatticeData assumes that the index for GetFOld is in distribution-space not site-space
           * @param distributionIndex
           * @return
           */
          distribn_t* GetFOld(site_t distributionIndex)
          {
            site_t globalIndex = distributionIndex / latticeInfo.GetNumVectors();
            site_t direction = distributionIndex % latticeInfo.GetNumVectors();
            return &distributions[globalIndex][direction];
          }
          /**
           * Get a pointer to the fOld array starting at the requested index. This version
           * of the function allows us to access the fOld array in a const way from a const
           * LatticeData
           * @param distributionIndex
           * @return
           */
          const distribn_t* GetFOld(site_t distributionIndex) const
          {
            site_t globalIndex = distributionIndex / latticeInfo.GetNumVectors();
            site_t direction = distributionIndex % latticeInfo.GetNumVectors();
            return &distributions.find(globalIndex)->second[direction];
          }

          /*
           * This is not defined for Neighbouring Data.
           * Data streamed across boundaries is handled by the existing mechanism.
           */
          template<typename LatticeType>
          site_t GetStreamedIndex(site_t globalIndex, unsigned int iDirectionIndex) const
          {
            return -1;
          }

          template<typename LatticeType>
          double GetCutDistance(site_t globalIndex, int direction) const
          {
            return distanceToWall.find(globalIndex)->second[direction - 1];
          }

          /**
           * Get the site data object for the given index.
           * @param iSiteIndex
           * @return
           */
          inline SiteData GetSiteData(site_t globalIndex) const
          {
            return siteData.find(globalIndex)->second;
          }

        private:
          std::map<site_t, std::vector<distribn_t> > distributions; //! The distribution values for the previous time step
          std::map<site_t, std::vector<distribn_t> > distanceToWall; //! Hold the distance to the wall for each fluid site and direction
          std::map<site_t, util::Vector3D<distribn_t> > wallNormalAtSite; //! Holds the wall normal near the fluid site, where appropriate
          std::map<site_t, SiteData> siteData; //! Holds the SiteData for each site.
          const lb::lattices::LatticeInfo& latticeInfo;
      };

    }
  }
}

#endif
