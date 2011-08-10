#ifndef HEMELB_VIS_VIEWPOINT_H
#define HEMELB_VIS_VIEWPOINT_H

#include "vis/Location.h"

namespace hemelb
{
  namespace vis
  {
    class Viewpoint
    {
    public:
      Viewpoint();
 
      Location<float> RotateToViewpoint(const Location<float>& iVector) const;

      Location <float> Project(const Location<float>& p1) const;

      void SetViewpointPosition(float longitude,
				float latitude,
				const Location<float>& iLocalCentre,
				float rad,
				float distance);

      const Location<float>& GetViewpointCentre() const;

    private:
      static Location <float> Rotate(
	float iSinX,
	float iCosX,
	float iSinY,
	float iCosY,
	const Location<float>& iVector);
	  
      float mDistance;
      float SinYRotation, CosYRotation;
      float SinXRotation, CosXRotation;
      
      Location<float> mViewpointCentre;
    };
  }
}

#endif /* HEMELB_VIS_VIEWPOINT_H */
