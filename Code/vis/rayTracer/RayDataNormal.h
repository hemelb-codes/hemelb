// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

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

          static MPI_Datatype GetMPIType()
          {
            const int lRayDataNormalCount = 14;
            int lRayDataNormalBlocklengths[lRayDataNormalCount] = { 1,
                                                                    1,
                                                                    1,
                                                                    1,
                                                                    1,
                                                                    1,
                                                                    1,
                                                                    1,
                                                                    1,
                                                                    1,
                                                                    1,
                                                                    1,
                                                                    1,
                                                                    1 };
            MPI_Datatype lRayDataNormalTypes[lRayDataNormalCount] = { MPI_LB,
                                                                      MpiDataType<int> (),
                                                                      MpiDataType<int> (),
                                                                      MpiDataType<float> (),
                                                                      MpiDataType<float> (),
                                                                      MpiDataType<float> (),
                                                                      MpiDataType<float> (),
                                                                      MpiDataType<float> (),
                                                                      MpiDataType<float> (),
                                                                      MpiDataType<float> (),
                                                                      MpiDataType<float> (),
                                                                      MpiDataType<float> (),
                                                                      MpiDataType<float> (),
                                                                      MPI_UB };

            MPI_Aint lRayDataNormalDisps[lRayDataNormalCount];

            RayDataNormal example[2];

            MPI_Get_address(&example[0], &lRayDataNormalDisps[0]);
            MPI_Get_address(&example[0].i, &lRayDataNormalDisps[1]);
            MPI_Get_address(&example[0].j, &lRayDataNormalDisps[2]);
            MPI_Get_address(&example[0].mLengthBeforeRayFirstCluster, &lRayDataNormalDisps[3]);
            MPI_Get_address(&example[0].mCumulativeLengthInFluid, &lRayDataNormalDisps[4]);
            MPI_Get_address(&example[0].mDensityAtNearestPoint, &lRayDataNormalDisps[5]);
            MPI_Get_address(&example[0].mStressAtNearestPoint, &lRayDataNormalDisps[6]);
            MPI_Get_address(&example[0].mVelR, &lRayDataNormalDisps[7]);
            MPI_Get_address(&example[0].mVelG, &lRayDataNormalDisps[8]);
            MPI_Get_address(&example[0].mVelB, &lRayDataNormalDisps[9]);
            MPI_Get_address(&example[0].mStressR, &lRayDataNormalDisps[10]);
            MPI_Get_address(&example[0].mStressG, &lRayDataNormalDisps[11]);
            MPI_Get_address(&example[0].mStressB, &lRayDataNormalDisps[12]);
            MPI_Get_address(&example[1], &lRayDataNormalDisps[13]);

            for (int index = lRayDataNormalCount - 1; index >= 0; index--)
            {
              lRayDataNormalDisps[index] -= lRayDataNormalDisps[0];
            }

            MPI_Datatype type;
            MPI_Type_struct(lRayDataNormalCount,
                            lRayDataNormalBlocklengths,
                            lRayDataNormalDisps,
                            lRayDataNormalTypes,
                            &type);
            return type;
          }

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

}

#endif // HEMELB_VIS_RAYTRACER_RAYDATANORMAL_H
