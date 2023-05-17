// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_STREAMERS_GUOZHENGSHIDELEGATE_H
#define HEMELB_LB_STREAMERS_GUOZHENGSHIDELEGATE_H

#include "lb/iolets/InOutLetVelocity.h"
#include "lb/streamers/BaseStreamerDelegate.h"
#include "geometry/neighbouring/RequiredSiteInformation.h"
#include "geometry/neighbouring/NeighbouringDataManager.h"
#include "util/Vector3D.h"

namespace hemelb::lb::streamers
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
        using CollisionType = CollisionImpl;
        using LatticeType = typename CollisionType::CKernel::LatticeType;
        using const_span = typename LatticeType::const_span;

        GuoZhengShiDelegate(CollisionType& delegatorCollider, InitParams& initParams) :
              collider(delegatorCollider),
                  neighbouringLatticeData(initParams.latDat->GetNeighbouringData()),
                  bValues(initParams.boundaryObject), bbDelegate(delegatorCollider, initParams)
        {
            // Want to loop over each site this streamer is responsible for,
            // as specified in the siteRanges.
            for (auto [first, last]: initParams.siteRanges)
            {
              for (site_t localIndex = first; localIndex < last; ++localIndex)
              {
                geometry::Site<const geometry::Domain> localSite =
                    initParams.latDat->GetSite(localIndex);

                const LatticeVector& localSiteLocation = localSite.GetGlobalSiteCoords();

                // Ignore ones that aren't walls;
                if (!localSite.IsWall())
                {
                  hemelb::log::Logger::Log<hemelb::log::Error, hemelb::log::OnePerCore>("GZS streamer initialised with non wall site [%d, %d, %d]",
                                                                                        localSiteLocation.x(),
                                                                                        localSiteLocation.y(),
                                                                                        localSiteLocation.z());
                  continue;
                }

                // Iterate over every neighbouring direction from here.
                for (Direction direction = 0; direction < LatticeType::NUMVECTORS; ++direction)
                {
                  // If there's a wall in direction i, then we need the site in opposite direction
                  // So if there's no wall, there's no init needed.
                  if (!localSite.HasWall(direction))
                    continue;

                  Direction opp = LatticeType::INVERSEDIRECTIONS[direction];

                  // If there's a wall or an iolet in this direction, we will have to do something else.
                  if (localSite.HasWall(opp) || localSite.HasIolet(opp))
                    continue;

                  // We will need this site's data - work out what task it's data is on.
                  const LatticeVector neighbourLocation = localSiteLocation
                      + LatticeVector(LatticeType::CX[opp],
                                      LatticeType::CY[opp],
                                      LatticeType::CZ[opp]);
                  proc_t neighbourSiteHomeProc =
                      initParams.latDat->GetProcIdFromGlobalCoords(neighbourLocation);

                  // A solid site - this should have been picked up above by HasWall/HasIolet
                  if (neighbourSiteHomeProc == SITE_OR_BLOCK_SOLID)
                  {
                    hemelb::log::Logger::Log<hemelb::log::Error, hemelb::log::OnePerCore>("Inconsistent cut links/neighbour status for site [%d, %d, %d]",
                                                                                          localSiteLocation.x(),
                                                                                          localSiteLocation.y(),
                                                                                          localSiteLocation.z());
                    continue;
                  }
                  // If it's on this task, we don't need to request its data.
                  if (neighbourSiteHomeProc == initParams.latDat->GetLocalRank())
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

          /*
           * The outline of this method is as follows:
           *
           * if (wallDistance < 0.75)
           *   if (site.HasIolet(i))
           *     if (is not velocity iolet)
           *       Do SBB
           *     else
           *       Do Modified GZS2
           *   else
           *     if (site.HasWall(i))
           *       Do SBB
           *     else
           *       Do Regular GZS2
           * else
           *   Do GZS1
           */
          inline void StreamLink(const LbmParameters* lbmParams,
                                 geometry::FieldData& latDat,
                                 const geometry::Site<geometry::FieldData>& site,
                                 HydroVars<typename CollisionType::CKernel>& hydroVars,
                                 const Direction& iPrime)
          {
            Direction i = LatticeType::INVERSEDIRECTIONS[iPrime];
            // Get the distance to the boundary.
            double wallDistance = site.GetWallDistance<LatticeType>(iPrime);

            // Set up for GZS - do the extrapolation from this site - u_w1

            // Now we work out the hypothetical velocity of the solid site on the other side
            // of the wall.
            // Assume that the wall velocity (0) is linearly interpolated along the line
            // between the nearest fluid site and the solid site inside the wall.
            // Then 0 = velocityWall * wallDistance + velocityFluid * (1 - wallDistance)
            // Hence velocityWall = velocityFluid * (1 - 1/wallDistance)
            FVector<LatticeType> fWall;
            HydroVars<typename CollisionType::CKernel> hydroVarsWall(fWall);

            hydroVarsWall.density = hydroVars.density;
            hydroVarsWall.tau = hydroVars.tau;
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

            // When this can't be done (i.e. when there's a wall/iolet in the way), fall back to SBB
            if (wallDistance < 0.75)
            {
              if (site.HasIolet(i))
              {
                int boundaryId = site.GetIoletId();
                auto iolet = dynamic_cast<InOutLetVelocity*>(bValues->GetLocalIolet(boundaryId));
                if (iolet == nullptr)
                {
                  // SBB
                  return bbDelegate.StreamLink(lbmParams, latDat, site, hydroVars, iPrime);
                }
                else
                {
                  // Modified GZS - there is a velocity iolet blocking the neighbouring
                  // site who's data we would use for the second extrapolation.
                  // Use the imposed condition instead.
                  LatticePosition sitePos(site.GetGlobalSiteCoords());

                  LatticePosition neighPos(sitePos);
                  neighPos += LatticeType::CD[i];

                  LatticeVelocity neighbourVelocity(iolet->GetVelocity(neighPos,
                                                                       bValues->GetTimeStep()));

                  // Obtain a second estimate, this time ignoring the fluid site closest to
                  // the wall. Interpolating the next site away and the site within the wall
                  // to the point on the wall itself (velocity 0):
                  // 0 = velocityWall * (1 + wallDistance) / 2 + velocityNextFluid * (1 - wallDistance)/2
                  // Rearranging gives velocityWall = velocityNextFluid * (wallDistance - 1)/(wallDistance+1)
                  LatticeVelocity velocityWallSecondEstimate = neighbourVelocity
                      * (wallDistance - 1) / (wallDistance + 1);
                  // Next, we interpolate between the first and second estimates to improve the estimate.
                  // Extrapolate to obtain the velocity at the wall site.
                  for (int dimension = 0; dimension < 3; dimension++)
                  {
                    hydroVarsWall.momentum[dimension] = wallDistance
                        * hydroVarsWall.momentum[dimension]
                        + (1. - wallDistance) * hydroVars.density
                            * velocityWallSecondEstimate[dimension];
                  }
                  // Should interpolate in the same way to get f_neq - skip since not available
                }
              }
              else
              {
                if (site.HasWall(i))
                {
                  // SBB
                  return bbDelegate.StreamLink(lbmParams, latDat, site, hydroVars, iPrime);
                }
                else
                {
                  // There is a neighbour site to use for standard GZS to calculate u_w2.
                  auto neighbourFOld = GetNeighbourFOld(site, i, latDat);
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
                                                             neighbourMomentum,
                                                             neighbourVelocity,
                                                             neighbourFEq);
                  }
                  // Obtain a second estimate, this time ignoring the fluid site closest to
                  // the wall. Interpolating the next site away and the site within the wall
                  // to the point on the wall itself (velocity 0):
                  // 0 = velocityWall * (1 + wallDistance) / 2 + velocityNextFluid * (1 - wallDistance)/2
                  // Rearranging gives velocityWall = velocityNextFluid * (wallDistance - 1)/(wallDistance+1)
                  LatticeVelocity velocityWallSecondEstimate = neighbourVelocity
                      * (wallDistance - 1) / (wallDistance + 1);
                  // Next, we interpolate between the first and second estimates to improve the estimate.
                  // Extrapolate to obtain the velocity at the wall site.
                  for (int dimension = 0; dimension < 3; dimension++)
                  {
                    hydroVarsWall.momentum[dimension] = wallDistance
                        * hydroVarsWall.momentum[dimension]
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
              }
            }
            else
            {
              // We don't need any neighbours so this works any way.
            }
            // Finally, we want to collide and stream, using the chosen collision kernel.
            //
            // hydroVarsWall contains the correct density and momentum, but
            // we must also set f, f_eq and f_neq to ensure that the collision
            // will work properly for more complex operators.

            // Calculate equilibrium values
            LatticeType::CalculateFeq(hydroVarsWall.density,
                                      hydroVarsWall.momentum,
                                      hydroVarsWall.GetFEq());

            // For the wall site, construct f_old  = f_eq + f_neq
            for (unsigned j = 0; j < LatticeType::NUMVECTORS; ++j)
            {
              fWall[j] = hydroVarsWall.GetFEq()[j] + hydroVarsWall.GetFNeq()[j];
            }
            // Perform collision
            collider.Collide(lbmParams, hydroVarsWall);
            // stream
            distribn_t* fNew = latDat.GetFNew(site.GetIndex() * LatticeType::NUMVECTORS);
            fNew[i] = hydroVarsWall.GetFPostCollision()[i];

          }

    private:
        const_span GetNeighbourFOld(const geometry::Site<geometry::FieldData>& site,
                                    const Direction& i,
                                    geometry::FieldData& latDat)
        {
            auto&& domain = latDat.GetDomain();
            // Find the neighbour's global location and which proc it's on.
            LatticeVector neighbourGlobalLocation = site.GetGlobalSiteCoords() + LatticeType::VECTORS[i];
            proc_t neighbourProcessor = domain.GetProcIdFromGlobalCoords(neighbourGlobalLocation);
            if (neighbourProcessor == domain.GetLocalRank())
            {
                // If it's local, get a Site object for it.
                geometry::Site<geometry::FieldData> nextSiteOut =
                        latDat.GetSite(domain.GetContiguousSiteId(neighbourGlobalLocation));
                return nextSiteOut.GetFOld<LatticeType>();
            }
            else
            {
                auto neighbourSite = latDat.GetNeighbouringData().GetSite(
                        domain.GetGlobalNoncontiguousSiteIdFromGlobalCoords(neighbourGlobalLocation)
                );
                return neighbourSite.GetFOld<LatticeType>();
            }
        }

        // the collision
        CollisionType collider;
        const geometry::neighbouring::NeighbouringDomain& neighbouringLatticeData;
        BoundaryValues* bValues;
        SimpleBounceBackDelegate<CollisionType> bbDelegate;
    };

}

#endif /* HEMELB_LB_STREAMERS_GUOZHENGSHIDELEGATE_H */
