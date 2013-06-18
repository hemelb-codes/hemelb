// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_LB_STREAMERS_GUOZHENGSHIDELEGATE_H
#define HEMELB_LB_STREAMERS_GUOZHENGSHIDELEGATE_H

#include "lb/streamers/BaseStreamerDelegate.h"
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
      class GuoZhengShiDelegate : public BaseStreamerDelegate<CollisionImpl>
      {
        public:
          typedef CollisionImpl CollisionType;
          typedef typename CollisionType::CKernel::LatticeType LatticeType;

          GuoZhengShiDelegate(CollisionType& delegatorCollider, kernels::InitParams& initParams) :
              collider(delegatorCollider),
                  neighbouringLatticeData(initParams.latDat->GetNeighbouringData())
          {
            // Want to loop over each site this streamer is responsible for,
            // as specified in the siteRanges.
            for (std::vector<std::pair<site_t, site_t> >::iterator rangeIt =
                initParams.siteRanges.begin(); rangeIt != initParams.siteRanges.end(); ++rangeIt)
            {
              for (site_t localIndex = rangeIt->first; localIndex < rangeIt->second; ++localIndex)
              {
                geometry::Site<const geometry::LatticeData> localSite =
                    initParams.latDat->GetSite(localIndex);

                // Ignore ones that aren't walls;
                if (!localSite.IsWall())
                  continue;

                const LatticeVector localSiteLocation = localSite.GetGlobalSiteCoords();

                // Iterate over every neighbouring direction from here.
                for (Direction direction = 0; direction < LatticeType::NUMVECTORS; ++direction)
                {
                  const LatticeVector neighbourLocation = localSiteLocation
                      + LatticeVector(LatticeType::CX[direction],
                                      LatticeType::CY[direction],
                                      LatticeType::CZ[direction]);

                  // Make sure we don't try to get info about off-lattice neighbours.
                  if (!initParams.latDat->IsValidLatticeSite(neighbourLocation))
                    continue;

                  proc_t neighbourSiteHomeProc =
                      initParams.latDat->GetProcIdFromGlobalCoords(neighbourLocation);

                  // BIG_NUMBER2 means a solid site. We don't want info about solids or
                  // neighbouring sites on this proc.
                  if (neighbourSiteHomeProc == BIG_NUMBER2
                      || neighbourSiteHomeProc
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
          }

          inline void StreamLink(const LbmParameters* lbmParams,
                                 geometry::LatticeData* const latDat,
                                 const geometry::Site<geometry::LatticeData>& site,
                                 kernels::HydroVars<typename CollisionType::CKernel>& hydroVars,
                                 const Direction& iPrime)
          {
            Direction i = LatticeType::INVERSEDIRECTIONS[iPrime];
            // Get the distance to the boundary.
            double wallDistance = site.GetWallDistance<LatticeType>(iPrime);
            // Now we work out the hypothetical velocity of the solid site on the other side
            // of the wall.
            // Assume that the wall velocity (0) is linearly interpolated along the line
            // between the nearest fluid site and the solid site inside the wall.
            // Then 0 = velocityWall * wallDistance + velocityFluid * (1 - wallDistance)
            // Hence velocityWall = velocityFluid * (1 - 1/wallDistance)
            distribn_t fWall[LatticeType::NUMVECTORS];
            kernels::HydroVars<typename CollisionType::CKernel> hydroVarsWall(fWall);

            hydroVarsWall.density = hydroVars.density;
            hydroVarsWall.momentum = hydroVars.momentum * (1. - 1. / wallDistance);

            // Find the non-equilibrium distribution in the unstreamed direction.
            std::copy(hydroVars.GetFNeqPtr(),
                      hydroVars.GetFNeqPtr() + LatticeType::NUMVECTORS,
                      hydroVarsWall.GetFNeqPtr());

            // The authors suggest simply using the wall velocity when there's a large distance
            // between the fluid site and wall (> 0.75 lattice vector).
            // When there's a smaller distance, they recommend looking at the next fluid site out
            // and extrapolating from that to obtain another estimate. The wallDistance is then used
            // as an interpolation variable between the two estimates.
            // A similar thing is done with the non-equilibrium distribution estimate. It is either
            // the value in that direction at the nearest site, or an interpolation between the values
            // at the nearest site and the next site away.
            if (wallDistance < 0.75 && !site.HasWall(i) && !site.HasIolet(i))
            {
              // We can only do this if there's gonna be a point there to interpolate from, i.e. there's no boundary
              // in direction i
              // Need some info about the next node away from the wall in this direction...
              // Populate this depending whether it's on this core or not.
              const distribn_t *neighbourFOld = GetNeighbourFOld(site, i, latDat);
              // Now calculate this field information.
              LatticeVelocity neighbourVelocity;
              distribn_t neighbourFEq[LatticeType::NUMVECTORS];
              // Go ahead and calculate the density, momentum and eqm distribution.
              {
                distribn_t neighbourDensity;
                LatticeVelocity neighbourMomentum;
                // Note that nextNodeOutVelocity is passed as the momentum argument, this
                // is because it is immediately divided by density when the function returns.
                LatticeType::CalculateDensityMomentumFEq(neighbourFOld,
                                                         neighbourDensity,
                                                         neighbourMomentum.x,
                                                         neighbourMomentum.y,
                                                         neighbourMomentum.z,
                                                         neighbourVelocity.x,
                                                         neighbourVelocity.y,
                                                         neighbourVelocity.z,
                                                         neighbourFEq);
              }
              // Obtain a second estimate, this time ignoring the fluid site closest to
              // the wall. Interpolating the next site away and the site within the wall
              // to the point on the wall itself (velocity 0):
              // 0 = velocityWall * (1 + wallDistance) / 2 + velocityNextFluid * (1 - wallDistance)/2
              // Rearranging gives velocityWall = velocityNextFluid * (wallDistance - 1)/(wallDistance+1)
              LatticeVelocity velocityWallSecondEstimate = neighbourVelocity * (wallDistance - 1)
                  / (wallDistance + 1);
              // Next, we interpolate between the first and second estimates to improve the estimate.
              // Extrapolate to obtain the velocity at the wall site.
              for (int dimension = 0; dimension < 3; dimension++)
              {
                hydroVarsWall.momentum[dimension] = wallDistance * hydroVarsWall.momentum[dimension]
                    + (1. - wallDistance) * hydroVars.density
                        * velocityWallSecondEstimate[dimension];
              }
              // Interpolate in the same way to get f_neq.
              distribn_t* fNeqWall = hydroVarsWall.GetFNeqPtr();
              for (unsigned j = 0; j < LatticeType::NUMVECTORS; ++j)
              {
                fNeqWall[j] = wallDistance * fNeqWall[j]
                    + (1. - wallDistance) * (neighbourFOld[j] - neighbourFEq[j]);
              }
            }

            // Finally, we want to collide and stream, using the chosen collision kernel.
            //
            // hydroVarsWall contains the correct density and momentum, but
            // we must also set f, f_eq and f_neq to ensure that the collision
            // will work properly for more complex operators.

            // Calculate equilibrium values
            LatticeType::CalculateFeq(hydroVarsWall.density,
                                      hydroVarsWall.momentum.x,
                                      hydroVarsWall.momentum.y,
                                      hydroVarsWall.momentum.z,
                                      hydroVarsWall.GetFEqPtr());

            // For the wall site, construct f_old  = f_eq + f_neq
            for (unsigned j = 0; j < LatticeType::NUMVECTORS; ++j)
            {
              fWall[j] = hydroVarsWall.GetFEq()[j] + hydroVarsWall.GetFNeq()[j];
            }
            // Perform collision
            collider.Collide(lbmParams, hydroVarsWall);
            // stream
            distribn_t* fNew = latDat->GetFNew(site.GetIndex() * LatticeType::NUMVECTORS);
            fNew[i] = hydroVarsWall.GetFPostCollision()[i];
          }

        private:
          const distribn_t *GetNeighbourFOld(const geometry::Site<geometry::LatticeData>& site,
                                             const Direction& i,
                                             geometry::LatticeData* const latDat)
          {
            const distribn_t* neighbourFOld;
            // Find the neighbour's global location and which proc it's on.
            LatticeVector neighbourGlobalLocation = site.GetGlobalSiteCoords()
                + LatticeVector(LatticeType::CX[i], LatticeType::CY[i], LatticeType::CZ[i]);
            proc_t neighbourProcessor = latDat->GetProcIdFromGlobalCoords(neighbourGlobalLocation);
            if (neighbourProcessor == topology::NetworkTopology::Instance()->GetLocalRank())
            {
              // If it's local, get a Site object for it.
              geometry::Site<geometry::LatticeData> nextSiteOut =
                  latDat->GetSite(latDat->GetContiguousSiteId(neighbourGlobalLocation));
              neighbourFOld = nextSiteOut.GetFOld<LatticeType>();
            }
            else
            {
              const geometry::neighbouring::ConstNeighbouringSite neighbourSite =
                  neighbouringLatticeData.GetSite(latDat->GetGlobalNoncontiguousSiteIdFromGlobalCoords(neighbourGlobalLocation));
              neighbourFOld = neighbourSite.GetFOld<LatticeType>();
            }
            return neighbourFOld;

          }
          // the collision
          CollisionType collider;
          const geometry::neighbouring::NeighbouringLatticeData& neighbouringLatticeData;
      };

    }
  }
}

#endif /* HEMELB_LB_STREAMERS_GUOZHENGSHIDELEGATE_H */
