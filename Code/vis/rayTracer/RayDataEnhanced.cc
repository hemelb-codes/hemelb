//#define NDEBUG;
#include <assert.h>
#include <iostream>

#include "vis/rayTracer/RayDataEnhanced.h"
#include "vis/rayTracer/HSLToRGBConverter.h"
#include "vis/ColPixel.h"
#include "vis/DomainStats.h"
#include "vis/Vector3D.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer
    {
      const float RayDataEnhanced::mMinLogIntensityMultiple = 0.5F;

      RayDataEnhanced::RayDataEnhanced() :
	mVelocitySum(0.0F),
	mStressSum(0.0F),
	mIntensity(0.0F)
      {
      }

      void RayDataEnhanced::DoUpdateData(const SiteData_t& iSiteData, 
					 const double* iWallNormal, 
					 const Vector3D<float>& iRayDirection,
					 const float iRayLengthInVoxel,
					 const float iAbsoluteDistanceFromViewpoint,
					 const DomainStats& iDomainStats,
					 const VisSettings& iVisSettings)
      {
	//assert(iSiteData.Density >= 0.0F);
	
	mVelocitySum += iSiteData.Velocity * (float) iDomainStats.velocity_threshold_max_inv;
	
	if (iVisSettings.mStressType != lb::ShearStress)
	{
	  // update the volume rendering of the von Mises stress flow field
	  mStressSum = iSiteData.Stress * (float) iDomainStats.stress_threshold_max_inv;
	}

        //Calculate the absolute dot product of the wall
	//vector normal and the ray direction
	float lDotProduct =
	  fabs(iRayDirection.DotProduct(
		 Vector3D<float>(iWallNormal[0],
				 iWallNormal[1],
				 iWallNormal[2])));

	//We want the attentuation to be between
	//mMinLogIntensityMultiple and 1
	float lIntensityMultiple = mMinLogIntensityMultiple +
	  (1.0F-mMinLogIntensityMultiple)*lDotProduct;

	//Scale the current intensity of the ray
	mIntensity *= lIntensityMultiple; 
      }
      
      void RayDataEnhanced::DoGetVelocityColour(unsigned char oColour[3]) const
      {
	HSLToRGBConverter::ConvertHSLToRGB(GetAverageVelocity()*180.0F,
					   1.0F,
					   GetIntensity(),
					   oColour);
      }
    

      void RayDataEnhanced::DoGetStressColour(unsigned char oColour[3]) const
      {
	HSLToRGBConverter::ConvertHSLToRGB(0.0F,
					   GetAverageStress(),
					   GetIntensity(),
					   oColour);
      }

     
      void RayDataEnhanced::DoMergeIn(const RayDataEnhanced& iOtherRayData,
				    const VisSettings& iVisSettings)
      {
	mVelocitySum += iOtherRayData.GetVelocitySum();

	if (iVisSettings.mStressType != lb::ShearStress)
	{
	  mStressSum += iOtherRayData.GetStressSum();;
	}
	
	mIntensity *= iOtherRayData.GetIntensity();
      }
    
      float RayDataEnhanced::GetVelocitySum() const
      {
	return mVelocitySum;
      }
      
      float RayDataEnhanced::GetStressSum() const
      {
	return mStressSum;
      }
      
      float RayDataEnhanced::GetIntensity() const
      {
	return mIntensity;
      }

      float RayDataEnhanced::GetAverageVelocity() const
      {
	return GetVelocitySum()/GetCumulativeLengthInFluid();
      }
      
      float RayDataEnhanced::GetAverageStress() const
      {
	return GetStressSum()/GetCumulativeLengthInFluid();
      }
    }
  }

  template<>
  MPI_Datatype MpiDataTypeTraits<hemelb::vis::raytracer::RayDataEnhanced>::RegisterMpiDataType()
  {
    int lRayDataEnhancedCount = 8;
    int lRayDataEnhancedBlocklengths[8] = { 1, 1, 1, 1, 1, 1, 1, 1};
    MPI_Datatype lRayDataEnhancedTypes[8] = { MPI_FLOAT,
					      MPI_FLOAT,
					      MPI_FLOAT,
					      MPI_FLOAT,
					      MPI_FLOAT,
					      MPI_FLOAT,
					      MPI_FLOAT,
					      MPI_UB };

    MPI_Aint lRayDataEnhancedDisps[8];

    lRayDataEnhancedDisps[0] = 0;

    for (int i = 1; i < lRayDataEnhancedCount; i++)
    {
      if (lRayDataEnhancedTypes[i - 1] == MPI_FLOAT)
      {
	lRayDataEnhancedDisps[i] = lRayDataEnhancedDisps[i - 1]
	  + (sizeof(float) * lRayDataEnhancedBlocklengths[i - 1]);
      }
      else if (lRayDataEnhancedTypes[i - 1] == MPI_INT)
      {
	lRayDataEnhancedDisps[i] = lRayDataEnhancedDisps[i - 1] + (sizeof(int) * lRayDataEnhancedBlocklengths[i - 1]);
      }
      else if (lRayDataEnhancedTypes[i - 1] == MPI_UNSIGNED)
      {
	lRayDataEnhancedDisps[i] = lRayDataEnhancedDisps[i - 1] + (sizeof(unsigned) * lRayDataEnhancedBlocklengths[i - 1]);
      }
    }
    MPI_Datatype type;
    MPI_Type_struct(lRayDataEnhancedCount,
		    lRayDataEnhancedBlocklengths,
		    lRayDataEnhancedDisps,
		    lRayDataEnhancedTypes,
		    &type);
    MPI_Type_commit(&type);
    return type;
  }
}
