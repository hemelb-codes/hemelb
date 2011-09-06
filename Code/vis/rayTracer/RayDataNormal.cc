//#define NDEBUG;
#include <assert.h>
#include <iostream>

#include "vis/rayTracer/RayDataNormal.h"
#include "vis/ColPixel.h"
#include "vis/DomainStats.h"
#include "vis/Vector3D.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer
    {
      
      RayDataNormal::RayDataNormal()
      {
	mVelR = 0.0F;
	mVelG = 0.0F;
	mVelB = 0.0F;

	mStressR = 0.0F;
	mStressG = 0.0F;
	mStressB = 0.0F;
      }

      void RayDataNormal::DoUpdateData(const SiteData_t& iSiteData, 
				       const double* iWallNormal, 
				       const float iRayLengthInVoxel,
				       const float iAbsoluteDistanceFromViewpoint,
				       const DomainStats& iDomainStats,
				       const VisSettings& iVisSettings)
      {

	assert(iSiteData.Density >= 0.0F);
	
	float lPalette[3];

	// update the volume rendering of the velocity flow field
	ColPixel<RayDataNormal>::PickColour(iSiteData.Velocity * (float) iDomainStats.velocity_threshold_max_inv,
					    lPalette);

	UpdateVelocityColour(iRayLengthInVoxel, lPalette);

	if (iVisSettings.mStressType != lb::ShearStress)
	{
	  // update the volume rendering of the von Mises stress flow field
	  float lScaledStress = iSiteData.Stress * (float) iDomainStats.stress_threshold_max_inv;

	  ColPixel<RayDataNormal>::PickColour(lScaledStress, lPalette);

	  UpdateStressColour(iRayLengthInVoxel, lPalette);
	}
      }
      
      void RayDataNormal::DoGetVelocityColour(unsigned char oColour[3]) const
      {
	oColour[0] = static_cast <unsigned char> ((mVelR*255.0F) / GetCumulativeLengthInFluid());
	oColour[1] = static_cast <unsigned char> ((mVelG*255.0F) / GetCumulativeLengthInFluid());
	oColour[2] = static_cast <unsigned char> ((mVelB*255.0F) / GetCumulativeLengthInFluid());
      }
    

      void RayDataNormal::DoGetStressColour(unsigned char oColour[3]) const
      {
	oColour[0] = static_cast<unsigned char> ((mStressR*255.0F) / GetCumulativeLengthInFluid());
	oColour[1] = static_cast<unsigned char> ((mStressG*255.0F) / GetCumulativeLengthInFluid());
	oColour[2] = static_cast<unsigned char> ((mStressB*255.0F) / GetCumulativeLengthInFluid());
      }

     
      void RayDataNormal::DoMergeIn(const RayDataNormal& iOtherRayData,
				    const VisSettings& iVisSettings)
      {
	mVelR += iOtherRayData.mVelR;
	mVelG += iOtherRayData.mVelG;
	mVelB += iOtherRayData.mVelB;

	if (iVisSettings.mStressType != lb::ShearStress)
	{
	  mStressR += iOtherRayData.mStressR;
	  mStressG += iOtherRayData.mStressG;
	  mStressB += iOtherRayData.mStressB;
	}
      }
    
      void RayDataNormal::UpdateVelocityColour(float iDt, const float iPalette[3])
      {
	mVelR += iDt * iPalette[0];
	mVelG += iDt * iPalette[1];
	mVelB += iDt * iPalette[2];
      }

      void RayDataNormal::UpdateStressColour(float iDt, const float iPalette[3])
      {
	mStressR += iDt * iPalette[0];
	mStressG += iDt * iPalette[1];
	mStressB += iDt * iPalette[2];
      }     
    }
  }

  template<>
  MPI_Datatype MpiDataTypeTraits<hemelb::vis::raytracer::RayDataNormal>::RegisterMpiDataType()
  {
    int lRayDataNormalCount = 11;
    int lRayDataNormalBlocklengths[11] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    MPI_Datatype lRayDataNormalTypes[11] = { MPI_FLOAT,
					     MPI_FLOAT,
					     MPI_FLOAT,
					     MPI_FLOAT,
					     MPI_FLOAT,
					     MPI_FLOAT,
					     MPI_FLOAT,
					     MPI_FLOAT,
					     MPI_FLOAT,
					     MPI_FLOAT,
					     MPI_UB };

    MPI_Aint lRayDataNormalDisps[11];

    lRayDataNormalDisps[0] = 0;

    for (int i = 1; i < lRayDataNormalCount; i++)
    {
      if (lRayDataNormalTypes[i - 1] == MPI_FLOAT)
      {
	lRayDataNormalDisps[i] = lRayDataNormalDisps[i - 1]
	  + (sizeof(float) * lRayDataNormalBlocklengths[i - 1]);
      }
      else if (lRayDataNormalTypes[i - 1] == MPI_INT)
      {
	lRayDataNormalDisps[i] = lRayDataNormalDisps[i - 1] + (sizeof(int) * lRayDataNormalBlocklengths[i - 1]);
      }
      else if (lRayDataNormalTypes[i - 1] == MPI_UNSIGNED)
      {
	lRayDataNormalDisps[i] = lRayDataNormalDisps[i - 1] + (sizeof(unsigned) * lRayDataNormalBlocklengths[i - 1]);
      }
    }
    MPI_Datatype type;
    MPI_Type_struct(lRayDataNormalCount,
		    lRayDataNormalBlocklengths,
		    lRayDataNormalDisps,
		    lRayDataNormalTypes,
		    &type);
    MPI_Type_commit(&type);
    return type;
  }
}
