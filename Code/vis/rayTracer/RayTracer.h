#ifndef HEMELB_VIS_RAYTRACER_H
#define HEMELB_VIS_RAYTRACER_H

#include <map>
#include <stack>
#include <vector>

#include "constants.h"
#include "geometry/LatticeData.h"
#include "topology/NetworkTopology.h"

#include "vis/DomainStats.h"
#include "vis/Screen.h"
#include "vis/Viewpoint.h"
#include "vis/VisSettings.h"
#include "vis/Vector3D.h"
#include "vis/XYCoordinates.h"
#include "vis/rayTracer/Cluster.h"
#include "vis/rayTracer/ClusterBuilder.h"
#include "vis/rayTracer/SiteData.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer 
    {
      class RayTracer
      {
      public:
	// Constructor and destructor do all the usual stuff.
	RayTracer(const geometry::LatticeData* iLatDat,
		  const DomainStats* iDomainStats,
		  Screen* iScreen,
		  Viewpoint* iViewpoint,
		  VisSettings* iVisSettings);
	~RayTracer();

	//Calls the cluster builder to build the clusters
	//Must be explicitly called after construction as the method is 
	//virtual 
	virtual void BuildClusters();

	// Method to update the voxel corresponding to site i with its
	// newly calculated density, velocity and stress.
	void UpdateClusterVoxel(site_t i,
				distribn_t density,
				distribn_t velocity,
				distribn_t stress);

	// Render the current state into an image.
	void Render();

      protected:
	ClusterBuilder* mClusterBuilder;
	const geometry::LatticeData* mLatDat;

      private:
	struct Ray
	{
	  Vector3D<float> Direction;
	  Vector3D <float> InverseDirection;
	  float Length;

	  float VelocityColour[3];
	  float StressColour[3];
	  float Stress;
	  float Density;
	  float MinT;
	};

	struct AABB
	{
	  float acc_1, acc_2, acc_3, acc_4, acc_5, acc_6;
	};

	void RenderCluster(const Cluster& iCluster);

	void UpdateSubImageExtentForCorner
	  (const Vector3D<float>& iCorner,
	   XYCoordinates<float>& ioSubImageLowerLeft,
	   XYCoordinates<float>& ioSubImageUpperRight);
	  
	void UpdateRayData(const SiteData_t* iSiteData,
			   float ray_t,
			   float ray_segment,
			   Ray* bCurrentRay);

	void TraverseVoxels(const Vector3D<float>& block_min,
			    const Vector3D<float>& block_x,
			    const SiteData_t* iSiteData,
			    float t,
			    Ray* bCurrentRay,
			    const Vector3D<bool>& xyz_is_1);

	void TraverseBlocks(const Cluster* cluster,
			    const Vector3D<bool>& xyz_Is_1,
			    const Vector3D<float>& ray_dx,
			    Ray *bCurrentRay);

	void AABBvsRay(const AABB* aabb,
		       const Vector3D<float>& inverseDirection,
		       const Vector3D<bool>& xyzComponentIsPositive,
		       float* t_near,
		       float* t_far);

	void UpdateColour(float dt, const float palette[3], float col[3]);
	
	const DomainStats* mDomainStats;
	Screen* mScreen;
	Viewpoint* mViewpoint;
	VisSettings* mVisSettings;

	const float mBlockSizeFloat;
	const float mBlockSizeInverse;
	const site_t block_size2, block_size3, block_size_1;
	const site_t blocksYz;
      };
    }
  }
}

#endif // HEMELB_VIS_RAYTRACER_H
