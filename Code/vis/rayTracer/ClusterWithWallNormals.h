#ifndef HEMELB_VIS_RAYTRACER_CLUSTERWITHWALLNORMALS_H
#define HEMELB_VIS_RAYTRACER_CLUSTERWITHWALLNORMALS_H

#include "vis/rayTracer/ClusterShared.h"
#include "vis/rayTracer/SiteData.h"


namespace hemelb
{
  namespace vis
  {
    namespace raytracer 
    {
      class ClusterWithWallNormals : public ClusterShared <ClusterWithWallNormals>
	{
	public:
	  ClusterWithWallNormals();
	  
	  void DoResizeVectors();

	  void DoResizeVectorsForBlock(site_t iBlockNumber, site_t iSize);

	  double const* DoGetWallData(site_t iBlockNumber, site_t iSiteNumber) const;
	 
	  void DoSetWallData(site_t iBlockNumber, site_t iSiteNumber, double* iData);
       
	  static bool DoNeedsWallNormals();

	private:
	  std::vector<std::vector<double*> > WallNormals;
	};

    }
  }
}

#endif // HEMELB_VIS_RAYTRACER_CLUSTERWITHWALLNORMALS_H
