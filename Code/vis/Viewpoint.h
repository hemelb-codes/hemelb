#ifndef HEMELB_VIS_VIEWPOINT_H
#define HEMELB_VIS_VIEWPOINT_H

#include "vis/Vector3D.h"

namespace hemelb
{
  namespace vis
  {
    class Viewpoint
    {
    public:
      Viewpoint();
 
      Vector3D<float> RotateToViewpoint(const Vector3D<float>& iVector) const;

      Vector3D <float> Project(const Vector3D<float>& p1) const;

      void SetViewpointPosition(float longitude,
				float latitude,
				const Vector3D<float>& iLocalCentre,
				float rad,
				float distance);

      const Vector3D<float>& GetViewpointCentre() const;

    private:
      static Vector3D <float> Rotate(
	float iSinX,
	float iCosX,
	float iSinY,
	float iCosY,
	const Vector3D<float>& iVector);
	  
      float mDistance;
      float SinYRotation, CosYRotation;
      float SinXRotation, CosXRotation;
      
      Vector3D<float> mViewpointCentre;
    };
  }
}

#endif /* HEMELB_VIS_VIEWPOINT_H */
