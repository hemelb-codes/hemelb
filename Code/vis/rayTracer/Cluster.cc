#include "vis/rayTracer/Cluster.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer
    {
      Cluster::Cluster() { }

      void Cluster::ResizeVectors()
      {
	SiteData.resize(blocksX*blocksY*blocksZ);
      }
       
    }
  }
}
