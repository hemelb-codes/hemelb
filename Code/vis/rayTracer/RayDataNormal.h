// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_VIS_RAYTRACER_RAYDATANORMAL_H
#define HEMELB_VIS_RAYTRACER_RAYDATANORMAL_H

#include "vis/DomainStats.h"
#include "vis/rayTracer/RayData.h"
#include "vis/VisSettings.h"
#include "util/Vector3D.h"
#include "lb/LbmParameters.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer
    {
      // NB functions prefixed Do should only be called by the base class
      class RayDataNormal : public RayData<RayDataNormal>
      {
        public:
          RayDataNormal(int i, int j);
          RayDataNormal();

          // Used to process the ray data for a normal (non-wall) fluid site
          void DoUpdateDataForNormalFluidSite(const SiteData_t& iSiteData,
                                              const util::Vector3D<float>& iRayDirection,
                                              const float iRayLengthInVoxel,
                                              const VisSettings& iVisSettings);

          // Used to process the ray data for wall site
          void DoUpdateDataForWallSite(const SiteData_t& iSiteData,
                                       const util::Vector3D<float>& iRayDirection,
                                       const float iRayLengthInVoxel,
                                       const VisSettings& iVisSettings,
                                       const util::Vector3D<double>* iWallNormal);

          // Carries out the merging of the ray data in this
          // inherited type, for different segments of the same ray
          void DoCombine(const RayDataNormal& iOtherRayData);

          //Obtains the colour representing the velocity ray trace
          void DoGetVelocityColour(unsigned char oColour[3],
                                   const float iNormalisedDistanceToFirstCluster,
                                   const DomainStats& iDomainStats) const;

          // Obtains the colour representing the stress ray trace
          void DoGetStressColour(unsigned char oColour[3],
                                 const float iNormalisedDistanceToFirstCluster,
                                 const DomainStats& iDomainStats) const;

          void DoProcessTangentingVessel();

          static MPI_Datatype GetMPIType();

          // We need this because RayDataNormal uses it for every voxel update
          static const DomainStats* mDomainStats;

        private:
          void UpdateVelocityColour(float iDt, const float iPalette[3]);

          void UpdateStressColour(float iDt, const float iPalette[3]);

          void MakeColourComponent(float value, unsigned char& colour) const;

          float mVelR, mVelG, mVelB;
          float mStressR, mStressG, mStressB;
      };
    }
  }

  namespace net
  {
    template<>
    MPI_Datatype MpiDataTypeTraits<hemelb::vis::raytracer::RayDataNormal>::RegisterMpiDataType();
  }
}

#endif // HEMELB_VIS_RAYTRACER_RAYDATANORMAL_H
