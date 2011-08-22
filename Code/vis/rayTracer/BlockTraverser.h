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
      //The class also contains a record of which blocks have been visited, which
      //is neccessary for the algoritm which uses this. No locations are automatically
      //marked visited, and methods have been created to assist with random access
      //of the lattice data as required by the algorithm
      class BlockTraverser : public VolumeTraverser
      {
      public:
	BlockTraverser(const geometry::LatticeData * iLatDat);
	~BlockTraverser();

	site_t CurrentBlockNumber();

	Vector3D<site_t> GetSiteCoordinatesOfLowestSiteInCurrentBlock();

	//Tranverses the block until the next unvisited block is reached.
	//Returns false if the end of the Volume is reached
	bool GoToNextUnvisitedBlock();

	geometry::LatticeData::BlockData * GetCurrentBlockData();

	geometry::LatticeData::BlockData * GetBlockDataForLocation(const Vector3D<site_t>& iLocation);

	site_t GetBlockSize();

	SiteTraverser GetSiteTraverserForCurrentBlock();

	SiteTraverser GetSiteTraverserForLocation(const Vector3D<site_t>& iLocation);

	virtual site_t GetXCount();

	virtual site_t GetYCount();

	virtual site_t GetZCount();

	bool IsValidLocation(Vector3D<site_t> block);

	bool IsCurrentBlockVisited();

	bool IsBlockVisited(site_t n);
	bool IsBlockVisited(Vector3D<site_t> n);

	void MarkCurrentBlockVisited();

	void MarkBlockVisited(site_t n);
	void MarkBlockVisited(Vector3D<site_t> location);

      private:
	bool GoToNextBlock();

	const geometry::LatticeData * mLatticeData;

	bool* mBlockVisited;
      }; 
    }
  }
}


#endif // HEMELB_VIS_BLOCKTRAVERSER_H
