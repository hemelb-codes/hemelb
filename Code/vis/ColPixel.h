#ifndef HEMELB_VIS_COLPIXEL_H
#define HEMELB_VIS_COLPIXEL_H

//#define NDEBUG
#include <math.h>
#include <assert.h>
#include <iostream>

#include "mpiInclude.h"
#include "util/utilityFunctions.h"
#include "vis/DomainStats.h"
#include "vis/rayTracer/RayData.h"
#include "vis/rayTracer/RayDataNormal.h"
#include "vis/VisSettings.h"
#include "lb/LbmParameters.h"
namespace hemelb
{
  namespace vis
  {
    template<typename RayDataType>
      class ColPixel
    {
    public:
     ColPixel()
    {
      SetStreakline(false);
      SetGlyph(false);
    }
      

    ColPixel(int iI, int iJ, bool iIsGlyph)
      {
	SetGlyph(iIsGlyph);
	SetStreakline(false);
	SetI(iI);
	SetJ(iJ);
      }

    ColPixel(int iI, int iJ, float particleVelocity, float particleZ, int particleInletId)
      {
	SetGlyph(false);
	SetStreakline(true);
	SetI(iI);
	SetJ(iJ);
	
	particle_vel = particleVelocity;
	particle_z = particleZ;
	mParticleInletId = particleInletId;
      }

    ColPixel(int iI, int iJ, const RayDataType& iRayData) :
       mRayData(iRayData)     
      { 
	SetGlyph(false);
	SetStreakline(false);
	SetI(iI);
	SetJ(iJ);
	
      }

      /**
       * Merge data from the first ColPixel argument into the second
       * ColPixel argument.
       */
      void MergeIn(const ColPixel<RayDataType> &iOtherPixel, const VisSettings& iVisSettings)
      {
	MergeRayTracingData(iOtherPixel, iVisSettings);

	// Now merge glyph data
	// TODO - do we really need all this? 
	// It may be cheaper to merge in all cases.
	if (iVisSettings.mStressType != lb::ShearStress
	    && (iVisSettings.mode == VisSettings::ISOSURFACES
		|| iVisSettings.mode == VisSettings::ISOSURFACESANDGLYPHS))
	{
	  if (iOtherPixel.IsGlyph())
	  {
	    SetGlyph(true);
	  }
	}
	else
	{
#ifndef NO_STREAKLINES
	  // merge streakline data
	  if (iOtherPixel.IsStreakline())
	  {
	    if (!this->IsStreakline())
	    {
	      particle_z = iOtherPixel.particle_z;
	      particle_vel = iOtherPixel.particle_vel;
	      mParticleInletId = iOtherPixel.mParticleInletId;

	      SetStreakline(true);
	    }
	    else
	    {
	      if (iOtherPixel.particle_z < particle_z)
	      {
		particle_z = iOtherPixel.particle_z;
		particle_vel = iOtherPixel.particle_vel;
		mParticleInletId = iOtherPixel.mParticleInletId;
	      }
	    }
	  }
#endif
	}
      }

      void MergeRayTracingData(const ColPixel<RayDataType> &iOtherPixel, const VisSettings& iVisSettings)
      {
	if (iOtherPixel.ContainsRayData())
	{
	  if (this->ContainsRayData())
	  {
	    mRayData.MergeIn(iOtherPixel.GetRayData(), iVisSettings);
	  }
	  else
	  {
	    //Only the 'from' merge-pixel is ray-tracing
	    mRayData = iOtherPixel.GetRayData();
	  }
	}
      }

      void rawWritePixel(int *pixel_index,
			 unsigned char rgb_data[12],
			 const DomainStats* iDomainStats,
			 const VisSettings* visSettings) const
      {
	const int bits_per_char = sizeof(char) * 8;
	*pixel_index = (GetI() << (2 * bits_per_char)) + GetJ();

	if (ContainsRayData())
	{
	  // store velocity volume rendering colour
	  mRayData.GetVelocityColour(&rgb_data[0]);

	  if (visSettings->mStressType != lb::ShearStress)
	  {
	    // store von Mises stress volume rendering colour
	    mRayData.GetStressColour(&rgb_data[3]);
	  }
	  else if (mRayData.GetNearestStress() < NO_VALUE_F)
	  {
	    float stress_col[3];
	    PickColour(mRayData.GetNearestStress(), stress_col);

	    // store wall shear stress colour
	    MakePixelColour(int(255.0F * stress_col[0]),
			    int(255.0F * stress_col[1]),
			    int(255.0F * stress_col[2]),
			    &rgb_data[3]);
	  }
	  else
	  {
	    rgb_data[3] = rgb_data[4] = rgb_data[5] = 0;
	  }
	} // if (isRt)
	else
	{
	  for (int ii = 0; ii < 6; ++ii)
	  {
	    rgb_data[ii] = 255;
	  }
	}

	if (visSettings->mStressType != lb::ShearStress
	    && visSettings->mode == VisSettings::ISOSURFACES)
	{
	  float density_col[3], stress_col[3];
	  PickColour(mRayData.GetNearestDensity(), density_col);
	  PickColour(mRayData.GetNearestStress(), stress_col);

	  // store wall pressure colour
	  MakePixelColour(int(255.0F * density_col[0]),
			  int(255.0F * density_col[1]),
			  int(255.0F * density_col[2]),
			  &rgb_data[6]);

	  // store von Mises stress colour
	  MakePixelColour(int(255.0F * stress_col[0]),
			  int(255.0F * stress_col[1]),
			  int(255.0F * stress_col[2]),
			  &rgb_data[9]);

	}
	else if (visSettings->mStressType != lb::ShearStress
		 && visSettings->mode == VisSettings::ISOSURFACESANDGLYPHS)
	{
	  float density_col[3], stress_col[3];
	  PickColour(mRayData.GetNearestDensity(), density_col);
	  PickColour(mRayData.GetNearestStress(), stress_col);

	  if (ContainsRayData())
	  {
	    if (!IsGlyph())
	    {
	      density_col[0] += 1.0F;
	      density_col[1] += 1.0F;
	      density_col[2] += 1.0F;

	      stress_col[0] += 1.0F;
	      stress_col[1] += 1.0F;
	      stress_col[2] += 1.0F;
	    }

	    // store wall pressure (+glyph) colour
	    MakePixelColour(int(127.5F * density_col[0]),
			    int(127.5F * density_col[1]),
			    int(127.5F * density_col[2]),
			    &rgb_data[6]);

	    // store von Mises stress (+glyph) colour
	    MakePixelColour(int(127.5F * stress_col[0]),
			    int(127.5F * stress_col[1]),
			    int(127.5F * stress_col[2]),
			    &rgb_data[9]);
	  }
	  else
	  {
	    for (int ii = 6; ii < 12; ++ii)
	    {
	      rgb_data[ii] = 0;
	    }
	  }

	}
	else if (IsStreakline())
	{
	  float scaled_vel = (float) (particle_vel * iDomainStats->velocity_threshold_max_inv);
	  float particle_col[3];
	  PickColour(scaled_vel, particle_col);

	  // store particle colour
	  MakePixelColour(int(255.0F * particle_col[0]),
			  int(255.0F * particle_col[1]),
			  int(255.0F * particle_col[2]),
			  &rgb_data[6]);

	  for (int ii = 9; ii < 12; ++ii)
	  {
	    rgb_data[ii] = rgb_data[ii - 3];
	  }
	}
	else
	{
	  // store pressure colour
	  rgb_data[6] = rgb_data[7] = rgb_data[8] =
	    (unsigned char) util::NumericalFunctions::enforceBounds(
	      int(127.5F * mRayData.GetNearestDensity()),
	      0,
	      127);
	  
	  // store shear stress or von Mises stress
	  if (mRayData.GetNearestStress() <  NO_VALUE_F) 
	  {
	    rgb_data[9] = rgb_data[10] = rgb_data[11] =
	      (unsigned char) util::NumericalFunctions::enforceBounds(
		int(127.5F * mRayData.GetNearestStress()),
		0,
		127);
	  }
	  else
	  {
	    rgb_data[9] = rgb_data[10] = rgb_data[11] = 0;
	  }
	}
      }

      float GetDensity() const
      {
	return mRayData.GetNearestDensity();
      }

      float GetStress() const
      {
	return mRayData.GetNearestStress();
      }

      bool ContainsRayData() const
      {
	return mRayData.ContainsRayData();
      }

      const RayDataType& GetRayData() const
      {
	return mRayData;
      }

    void SetI(unsigned int iI)
      {
      	mI = iI;
      }

      void SetJ(unsigned int iJ)
      {
	mJ = iJ;
      }	  
	  
      unsigned int GetI() const
      {
	return mI;
      }

      unsigned int GetJ() const
      {
	return mJ;
      }

      void SetGlyph(bool iIsGlyph)
      {
	mIsGlyph = iIsGlyph;
      }

      unsigned int mIsGlyph;

      bool IsGlyph() const
      {
	return mIsGlyph;
      }

      void SetStreakline(bool iIsStreakline)
      {
	mIsStreakline = iIsStreakline;
      }
      unsigned int mIsStreakline;

      bool IsStreakline() const
      {
	return mIsStreakline;;
      }

/*
      void SetI(unsigned int iI)
      {
	assert(!(iI & mMostSignificantBit));
	mI = iI | (mI & mMostSignificantBit);
      }

      void SetJ(unsigned int iJ)
      {
	assert(!(iJ & mMostSignificantBit));
	mJ = iJ | (mJ && mMostSignificantBit);
      }	  
	  
      unsigned int GetI() const
      {
	return mI & ~mMostSignificantBit;
      }

      unsigned int GetJ() const
      {
	return mJ & ~mMostSignificantBit;
      }

      void SetGlyph(bool iIsGlyph)
      {
	if (iIsGlyph)
	{
	  mI = mI & mMostSignificantBit;
	}
	else
	{
	  mI = mI & ~mMostSignificantBit;
	}
      }


      bool IsGlyph() const
      {
	return mI && mMostSignificantBit;
      }

      void SetStreakline(bool iIsStreakline)
      {
	if (iIsStreakline)
	{
	  mJ = mJ & mMostSignificantBit;
	}
	else
	{
	  mJ = mJ & ~mMostSignificantBit;
	}
      }


      bool IsStreakline() const
      {
	return mJ && mMostSignificantBit;
      }*/

      static void PickColour(float value, float colour[3])
      {
	colour[0] = util::NumericalFunctions::enforceBounds<float>(4.F * value - 2.F, 0.F, 1.F);
	colour[1] = util::NumericalFunctions::enforceBounds<float>(2.F
								   - 4.F
								   * (float) fabs(value
										  - 0.5F),
								   0.F,
								   1.F);
	colour[2] = util::NumericalFunctions::enforceBounds<float>(2.F - 4.F * value, 0.F, 1.F);
      }

      static void MakePixelColour(int rawRed, int rawGreen, int rawBlue, unsigned char* dest)
      {
	dest[0] = (unsigned char) util::NumericalFunctions::enforceBounds(rawRed, 0, 255);
	dest[1] = (unsigned char) util::NumericalFunctions::enforceBounds(rawGreen, 0, 255);
	dest[2] = (unsigned char) util::NumericalFunctions::enforceBounds(rawBlue, 0, 255);
      }

    private:
      unsigned int mI; // Most significant bit used to indicate glyph
      unsigned int mJ; // Most significant bit used to indicated streakline

      static const unsigned int mMostSignificantBit = 1<<31;
	
      //Raytracing pixel data;
      RayDataType mRayData;

      // Streakline pixel data
      float particle_vel;
      float particle_z;
      int mParticleInletId;

      //NB please ensure that the MPI data type in
      //ColPixelMPIDataTypes.cc is maintained with
      //the private members
    };
  }
}

#endif // HEMELB_VIS_COLPIXEL_H
