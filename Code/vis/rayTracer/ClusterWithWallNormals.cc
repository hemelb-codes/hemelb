#include "geometry/LatticeData.h"
#include "vis/rayTracer/ClusterWithWallNormals.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer
    {
      ClusterWithWallNormals::ClusterWithWallNormals()
      { }

      void ClusterWithWallNormals::ResizeVectors()
      {
	SiteData.resize(blocksX*blocksY*blocksZ);
	WallNormals.resize(blocksX*blocksY*blocksZ);
      }
       
    }
  }
}
