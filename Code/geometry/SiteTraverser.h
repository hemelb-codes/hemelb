#ifndef HEMELB_GEOMETRY_SITETRAVERSER_H
#define HEMELB_GEOMETRY_SITETRAVERSER_H

#include "geometry/LatticeData.h"
#include "geometry/VolumeTraverser.h"

namespace hemelb
{
  namespace geometry
  {
    /**
     * SiteTraverser is used to traverse the sites in a speficied block
     * within the lattice data
     */
    class SiteTraverser : public VolumeTraverser
    {
      public:
        SiteTraverser(const geometry::LatticeData & iLatticeDat);

        virtual ~SiteTraverser();

        virtual site_t GetXCount();

        virtual site_t GetYCount();

        virtual site_t GetZCount();

      private:
        //Returns the block size in number of sites
        site_t GetBlockSize();

        const geometry::LatticeData & mLatticeData;
    };

  }
}

#endif // HEMELB_GEOMETRY_SITETRAVERSER_H
