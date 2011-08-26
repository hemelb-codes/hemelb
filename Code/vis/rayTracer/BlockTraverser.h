#ifndef HEMELB_VIS_BLOCKTRAVERSER_H
#define HEMELB_VIS_BLOCKTRAVERSER_H

#include "vis/rayTracer/VolumeTraverser.h"
#include "vis/rayTracer/SiteTraverser.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer 
    {
      //BlockTraverser is used to traverse the blocks in a lattice sequentially.
      class BlockTraverser : public VolumeTraverser
      {
      public:
	BlockTraverser(const geometry::LatticeData * iLatDat);
	~BlockTraverser();

	site_t CurrentBlockNumber();

	Vector3D<site_t> GetSiteCoordinatesOfLowestSiteInCurrentBlock();


	geometry::LatticeData::BlockData * GetCurrentBlockData();

	geometry::LatticeData::BlockData * GetBlockDataForLocation(const Vector3D<site_t>& iLocation);

	site_t GetBlockSize();

	SiteTraverser GetSiteTraverser();

	virtual site_t GetXCount();

	virtual site_t GetYCount();

	virtual site_t GetZCount();

	bool IsValidLocation(Vector3D<site_t> block);


      protected:
	bool GoToNextBlock();

	const geometry::LatticeData * mLatticeData;
      }; 
    }
  }
}


#endif // HEMELB_VIS_BLOCKTRAVERSER_H
