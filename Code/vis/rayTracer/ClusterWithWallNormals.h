#ifndef HEMELB_VIS_RAYTRACER_CLUSTERWITHWALLNORMALS_H
#define HEMELB_VIS_RAYTRACER_CLUSTERWITHWALLNORMALS_H

#include "vis/rayTracer/Cluster.h"
#include "vis/rayTracer/SiteData.h"


namespace hemelb
{
  namespace vis
  {
    namespace raytracer 
    {
      class ClusterWithWallNormals : public Cluster
	{
	public:
	  ClusterWithWallNormals();
	  ~ClusterWithWallNormals();

	  virtual void ResizeVectors();

	  std::vector<std::vector<double*> > WallNormals;
	};

    }
  }
}

#endif // HEMELB_VIS_RAYTRACER_CLUSTERWITHWALLNORMALS_H
