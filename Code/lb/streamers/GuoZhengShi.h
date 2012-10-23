// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_LB_STREAMERS_GUOZHENGSHI_H
#define HEMELB_LB_STREAMERS_GUOZHENGSHI_H

#include "lb/streamers/BaseStreamer.h"
#include "geometry/neighbouring/RequiredSiteInformation.h"
#include "geometry/neighbouring/NeighbouringDataManager.h"

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {
      /**
       * This class implements the boundary condition described by Guo, Zheng and Shi
       * in 'An Extrapolation Method for Boundary Conditions in Lattice-Boltzmann method'
       * Physics of Fluids, 14/6, June 2002, pp 2007-2010.
       */
      template<typename CollisionImpl>
      class GuoZhengShi : public BaseStreamer<GuoZhengShi<CollisionImpl> >
      {
        public:
          typedef CollisionImpl CollisionType;

        private:
          CollisionType collider;
          typedef typename CollisionType::CKernel::LatticeType LatticeType;

        public:
          GuoZhengShi(kernels::InitParams& initParams) :
            collider(initParams), neighbouringLatticeData(initParams.latDat->GetNeighbouringData())
          {
            // Go through every site on the local processor.
            for (site_t localIndex = 0; localIndex < initParams.latDat->GetLocalFluidSiteCount(); ++localIndex)
            {
              geometry::ConstSite localSite = initParams.latDat->GetSite(localIndex);

              // Ignore ones that aren't edges;
              if (!localSite.IsEdge())
                continue;

              const LatticeVector localSiteLocation = localSite.GetGlobalSiteCoords();

              // Iterate over every neighbouring direction from here.
              for (Direction direction = 0; direction < LatticeType::NUMVECTORS; ++direction)
              {
                const LatticeVector neighbourLocation = localSiteLocation + LatticeVector(LatticeType::CX[direction],
                                                                                          LatticeType::CY[direction],
                                                                                          LatticeType::CZ[direction]);

                // Make sure we don't try to get info about off-lattice neighbours.
                if (!initParams.latDat->IsValidLatticeSite(neighbourLocation))
                  continue;

                proc_t neighbourSiteHomeProc = initParams.latDat->GetProcIdFromGlobalCoords(neighbourLocation);

                // BIG_NUMBER2 means a solid site. We don't want info about solids or
                // neighbouring sites on this proc.
                if (neighbourSiteHomeProc == BIG_NUMBER2 || neighbourSiteHomeProc
                    == topology::NetworkTopology::Instance()->GetLocalRank())
                  continue;

                // Create a requirements with the info we need.
                geometry::neighbouring::RequiredSiteInformation requirements(false);

                requirements.Require(geometry::neighbouring::terms::Density);
                requirements.Require(geometry::neighbouring::terms::Velocity);

                initParams.neighbouringDataManager->RegisterNeededSite(initParams.latDat->GetGlobalNoncontiguousSiteIdFromGlobalCoords(neighbourLocation),
                                                                       requirements);
              }
            }
          }

          template<bool tDoRayTracing>
          inline void DoStreamAndCollide(const site_t firstIndex,
                                         const site_t siteCount,
                                         const LbmParameters* lbmParams,
                                         geometry::LatticeData* latDat,
                                         lb::MacroscopicPropertyCache& propertyCache)
          {
            for (site_t siteIndex = firstIndex; siteIndex < (firstIndex + siteCount); siteIndex++)
            {
              geometry::Site site = latDat->GetSite(siteIndex);

              // First do a normal collision & streaming step, as if we were mid-fluid.
              kernels::HydroVars<typename CollisionType::CKernel> hydroVars(site.GetFOld<LatticeType>());

              collider.CalculatePreCollision(hydroVars, site);
              collider.Collide(lbmParams, hydroVars);

              // Perform the streaming of the post-collision distribution.
              for (Direction direction = 0; direction < LatticeType::NUMVECTORS; direction++)
              {
                * (latDat->GetFNew(site.GetStreamedIndex<LatticeType>(direction))) =
                    hydroVars.GetFPostCollision()[direction];
              }

              // Now fill in the distributions that won't be streamed to
              // (those that point away from boundaries).
              // It's more convenient to iterate over the opposite direction.
              for (Direction oppositeDirection = 1; oppositeDirection < LatticeType::NUMVECTORS; oppositeDirection++)
              {
                // If there's a boundary in the opposite direction, we need to fill in the distribution.
                if (site.HasBoundary(oppositeDirection))
                {
                  Direction unstreamedDirection = LatticeType::INVERSEDIRECTIONS[oppositeDirection];

                  // Get the distance to the boundary.
                  double wallDistance = site.GetWallDistance<LatticeType> (oppositeDirection);

                  // Now we work out the hypothetical velocity of the solid site on the other side
                  // of the wall.
                  // Assume that the wall velocity (0) is linearly interpolated along the line
                  // between the nearest fluid site and the solid site inside the wall.
                  // Then 0 = velocityWall * wallDistance + velocityFluid * (1 - wallDistance)
                  // Hence velocityWall = velocityFluid * (1 - 1/wallDistance)
                  LatticeVelocity velocityWall = hydroVars.momentum * (1. - 1. / wallDistance) / hydroVars.density;

                  // Find the non-equilibrium distribution in the unstreamed direction.
                  distribn_t fNeqInUnstreamedDirection = hydroVars.GetFNeq()[unstreamedDirection];

                  // The authors suggest simply using the wall velocity when there's a large distance
                  // between the fluid site and wall (> 0.75 lattice vector).
                  // When there's a smaller distance, they recommend looking at the next fluid site out
                  // and extrapolating from that to obtain another estimate. The wallDistance is then used
                  // as an interpolation variable between the two estimates.
                  // A similar thing is done with the non-equilibrium distribution estimate. It is either
                  // the value in that direction at the nearest site, or an interpolation between the values
                  // at the nearest site and the next site away.
                  if (wallDistance < 0.75 && !site.HasBoundary(unstreamedDirection))
                  {
                    // We can only do this if there's gonna be a point there to interpolate from, i.e. there's no boundary
                    // in the direction of awayFromWallIndex.
                    // Need some info about the next node away from the wall in this direction...

                    // Find the neighbour's global location and which proc it's on.
                    LatticeVector neighbourGlobalLocation = site.GetGlobalSiteCoords()
                        + LatticeVector(LatticeType::CX[unstreamedDirection],
                                        LatticeType::CY[unstreamedDirection],
                                        LatticeType::CZ[unstreamedDirection]);

                    proc_t neighbourProcessor = latDat->GetProcIdFromGlobalCoords(neighbourGlobalLocation);

                    // Populate this depending whether it's on this core or not.
                    const distribn_t* neighbourFOld;

                    if (neighbourProcessor == topology::NetworkTopology::Instance()->GetLocalRank())
                    {
                      // If it's local, get a Site object for it.
                      geometry::Site nextSiteOut =
                          latDat->GetSite(latDat->GetContiguousSiteId(neighbourGlobalLocation));

                      neighbourFOld = nextSiteOut.GetFOld<LatticeType> ();
                    }
                    else
                    {
                      const geometry::neighbouring::ConstNeighbouringSite
                          neighbourSite =
                              neighbouringLatticeData.GetSite(latDat->GetGlobalNoncontiguousSiteIdFromGlobalCoords(neighbourGlobalLocation));

                      neighbourFOld = neighbourSite.GetFOld<LatticeType> ();
                    }

                    // Now calculate this field information.
                    LatticeVelocity nextNodeOutVelocity;
                    distribn_t nextNodeOutFEq[LatticeType::NUMVECTORS];

                    // Go ahead and calculate the density, momentum and eqm distribution.
                    {
                      distribn_t nextNodeDensity;
                      // Note that nextNodeOutVelocity is passed as the momentum argument, this
                      // is because it is immediately divided by density when the function returns.
                      LatticeType::CalculateDensityMomentumFEq(neighbourFOld,
                                                               nextNodeDensity,
                                                               nextNodeOutVelocity.x,
                                                               nextNodeOutVelocity.y,
                                                               nextNodeOutVelocity.z,
                                                               nextNodeOutFEq);
                      nextNodeOutVelocity /= nextNodeDensity;
                    }

                    // Obtain a second estimate, this time ignoring the fluid site closest to
                    // the wall. Interpolating the next site away and the site within the wall
                    // to the point on the wall itself (velocity 0):
                    // 0 = velocityWall * (1 + wallDistance) / 2 + velocityNextFluid * (1 - wallDistance)/2
                    // Rearranging gives velocityWall = velocityNextFluid * (wallDistance - 1)/(wallDistance+1)
                    LatticeVelocity velocityWallSecondEstimate = nextNodeOutVelocity * (wallDistance - 1)
                        / (wallDistance + 1);

                    // Next, we interpolate between the first and second estimates to improve the estimate.
                    // Extrapolate to obtain the velocity at the wall site.
                    for (int dimension = 0; dimension < 3; dimension++)
                    {
                      velocityWall[dimension] = wallDistance * velocityWall[dimension] + (1. - wallDistance)
                          * velocityWallSecondEstimate[dimension];
                    }

                    // Interpolate in the same way to get f_neq.
                    fNeqInUnstreamedDirection = wallDistance * fNeqInUnstreamedDirection + (1. - wallDistance)
                        * (neighbourFOld[unstreamedDirection] - nextNodeOutFEq[unstreamedDirection]);
                  }

                  // Use a helper function to calculate the actual value of f_eq in the desired direction at the wall node.
                  // Note that we assume that the density is the same as at this node
                  LatticeVelocity momentumWall = velocityWall * hydroVars.density;

                  distribn_t fEqTemp[LatticeType::NUMVECTORS];
                  LatticeType::CalculateFeq(hydroVars.density,
                                            momentumWall[0],
                                            momentumWall[1],
                                            momentumWall[2],
                                            fEqTemp);

                  // Collide and stream!
                  // TODO: It's not clear whether we should defer to the template collision type here
                  // or do a standard LBGK (implemented).
                  * (latDat->GetFNew(siteIndex * LatticeType::NUMVECTORS + unstreamedDirection))
                      = fEqTemp[unstreamedDirection] + (1.0 + lbmParams->GetOmega()) * fNeqInUnstreamedDirection;
                }
              }

              hydroVars.tau = lbmParams->GetTau();

              BaseStreamer<GuoZhengShi>::template UpdateMinsAndMaxes<tDoRayTracing>(site,
                                                                                    hydroVars,
                                                                                    lbmParams,
                                                                                    propertyCache);
            }
          }

          template<bool tDoRayTracing>
          inline void DoPostStep(const site_t firstIndex,
                                 const site_t siteCount,
                                 const LbmParameters* lbmParams,
                                 geometry::LatticeData* latDat,
                                 lb::MacroscopicPropertyCache& propertyCache)
          {

          }

        private:
          const geometry::neighbouring::NeighbouringLatticeData& neighbouringLatticeData;
      };
    }
  }
}

#endif /* HEMELB_LB_STREAMERS_GUOZHENGSHI_H */
