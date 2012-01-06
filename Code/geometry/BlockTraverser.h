#ifndef HEMELB_GEOMETRY_BLOCKTRAVERSER_H
#define HEMELB_GEOMETRY_BLOCKTRAVERSER_H

#include "lb/lattices/D3Q27.h"
#include "geometry/VolumeTraverser.h"
#include "geometry/SiteTraverser.h"

namespace hemelb
{
  namespace geometry
  {

    //TODO Ideally we'd hide implementation details like this.
    // Data about an element of the domain wall
    struct WallData
    {
        // estimated wall normal (if the site is close to the wall);
        double wall_nor[3];
        // cut distances along the 14 non-zero lattice vectors;
        // each one is between 0 and 1 if the surface cuts the corresponding
        // vector or is equal to "NO_VALUE" otherwise
        double cut_dist[lb::lattices::D3Q27::NUMVECTORS - 1];
    };

    // Data about each global block in the lattice,
    // site_data[] is an array containing individual lattice site data
    // within a global block.
    struct BlockData
    {
        BlockData()
        {
          ProcessorRankForEachBlockSite = NULL;
          wall_data = NULL;
          site_data = NULL;
        }

        ~BlockData()
        {
          if (ProcessorRankForEachBlockSite != NULL)
          {
            delete[] ProcessorRankForEachBlockSite;
            ProcessorRankForEachBlockSite = NULL;
          }
          if (wall_data != NULL)
          {
            delete[] wall_data;
            wall_data = NULL;
          }
          if (site_data != NULL)
          {
            delete[] site_data;
            site_data = NULL;
          }
        }

        // An array of the ranks on which each lattice site within the block resides.
        proc_t* ProcessorRankForEachBlockSite;
        // Information about wall / inlet / outlet position and orientation for
        // each site.
        WallData *wall_data;
        // The "site data" for each site.
        sitedata_t* site_data;
    };

    /**
     *BlockTraverser is used to traverse the blocks in a lattice sequentially.
     */
    class BlockTraverser : public VolumeTraverser
    {
      public:
        BlockTraverser(const LatticeData& iLatDat);
        ~BlockTraverser();

        site_t CurrentBlockNumber() const;

        util::Vector3D<site_t> GetSiteCoordinatesOfLowestSiteInCurrentBlock();

        BlockData* GetCurrentBlockData();

        BlockData* GetBlockDataForLocation(const util::Vector3D<site_t>& iLocation);

        site_t GetBlockSize();

        SiteTraverser GetSiteTraverser();

        site_t GetXCount() const;

        site_t GetYCount() const;

        site_t GetZCount() const;

        bool IsValidLocation(util::Vector3D<site_t> block);

      protected:
        bool GoToNextBlock();

        const LatticeData & mLatticeData;
    };
  }
}

#endif // HEMELB_GEOMETRY_BLOCKTRAVERSER_H
