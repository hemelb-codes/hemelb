
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_GEOMETRY_NEIGHBOURING_NEIGHBOURINGLATTICEDATA_H
#define HEMELB_GEOMETRY_NEIGHBOURING_NEIGHBOURINGLATTICEDATA_H
#include <map>
#include "geometry/Site.h"
#include "geometry/SiteData.h"
#include "lb/lattices/LatticeInfo.h"
namespace hemelb
{
  namespace geometry
  {
    namespace neighbouring
    {
      class NeighbouringSite;
      class ConstNeighbouringSite;

      // Here, all site indices are GLOBAL index.
      // Local users must determine the global index of the site they are interested in.

      class NeighbouringLatticeData
      {

        public:
          friend class Site<NeighbouringLatticeData> ; //! Let the inner classes have access to site-related data that's otherwise private.

          NeighbouringLatticeData(const lb::lattices::LatticeInfo& latticeInfo);
          virtual ~NeighbouringLatticeData()
          {
          }

          void SaveSite(site_t index,
                        const std::vector<distribn_t> &distribution,
                        const std::vector<distribn_t> &distances,
                        const util::Vector3D<distribn_t> &normal,
                        const SiteData & data);

          /**
           * Get a site object for the given index.
           * @param localIndex
           * @return
           */
          NeighbouringSite GetSite(site_t globalIndex);
          ConstNeighbouringSite GetSite(site_t globalIndex) const;
          /**
           * Get the wall normal at the given site
           * @param iSiteIndex
           * @return
           */
          const util::Vector3D<distribn_t>& GetNormalToWall(site_t globalIndex) const;
          util::Vector3D<distribn_t>& GetNormalToWall(site_t globalIndex);
          /**
           * Get a pointer to the fOld array starting at the requested index
           * LatticeData assumes that the index for GetFOld is in distribution-space not site-space
           * @param distributionIndex
           * @return
           */
          distribn_t* GetFOld(site_t distributionIndex);

          /**
           * Get a vector of the fOld array for the site
           * LatticeData assumes that the index for GetFOld is in distribution-space not site-space
           * @param globalIndex
           * @return
           */
          std::vector<distribn_t>& GetDistribution(site_t globalIndex);
          /**
           * Get a pointer to the fOld array starting at the requested index. This version
           * of the function allows us to access the fOld array in a const way from a const
           * LatticeData
           * @param distributionIndex
           * @return
           */
          const distribn_t* GetFOld(site_t distributionIndex) const;

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

          /*
           * For compatibility with lattice data,
           * these have to be distribn_t *, not a vector
           * because the lattice data stores the distances as a contiguous array
           */
          const distribn_t * GetCutDistances(site_t globalIndex) const;
          distribn_t* GetCutDistances(site_t globalIndex);

          /**
           * Get the site data object for the given index.
           * @param iSiteIndex
           * @return
           */
          const SiteData &GetSiteData(site_t globalIndex) const;
          SiteData &GetSiteData(site_t globalIndex);

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
