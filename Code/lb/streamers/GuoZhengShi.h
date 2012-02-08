#ifndef HEMELB_LB_STREAMERS_GUOZHENGSHI_H
#define HEMELB_LB_STREAMERS_GUOZHENGSHI_H

#include "lb/streamers/BaseStreamer.h"

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
      template<typename CollisionType>
      class GuoZhengShi : public BaseStreamer<GuoZhengShi<CollisionType> >
      {
        private:
          CollisionType collider;

        public:
          GuoZhengShi(kernels::InitParams& initParams) :
              collider(initParams)
          {

          }

          template<bool tDoRayTracing>
          inline void DoStreamAndCollide(const site_t firstIndex,
                                         const site_t siteCount,
                                         const LbmParameters* lbmParams,
                                         geometry::LatticeData* latDat,
                                         hemelb::vis::Control *visControl)
          {
            for (site_t siteIndex = firstIndex; siteIndex < (firstIndex + siteCount); siteIndex++)
            {
              geometry::Site site = latDat->GetSite(siteIndex);

              // First do a normal collision & streaming step, as if we were mid-fluid.
              kernels::HydroVars<typename CollisionType::CKernel> hydroVars(site.GetFOld());

              collider.CalculatePreCollision(hydroVars, site);
              collider.Collide(lbmParams, hydroVars);

              // Perform the streaming of the post-collision distribution.
              for (Direction direction = 0; direction < D3Q15::NUMVECTORS; direction++)
              {
                * (latDat->GetFNew(site.GetStreamedIndex(direction))) =
                    hydroVars.GetFPostCollision()[direction];
              }

              // Now fill in the distributions that won't be streamed to
              // (those that point away from boundaries).
              // It's more convenient to iterate over the opposite direction.
              for (Direction oppositeDirection = 1; oppositeDirection < D3Q15::NUMVECTORS;
                  oppositeDirection++)
              {
                // If there's a boundary in the opposite direction, we need to fill in the distribution.
                if (site.HasBoundary(oppositeDirection))
                {
                  Direction unstreamedDirection = D3Q15::INVERSEDIRECTIONS[oppositeDirection];

                  // Get the distance to the boundary.
                  double wallDistance = site.GetWallDistance(oppositeDirection);

                  // Now we work out the hypothetical velocity of the solid site on the other side
                  // of the wall.
                  // Assume that the wall velocity (0) is linearly interpolated along the line
                  // between the nearest fluid site and the solid site inside the wall.
                  // Then 0 = velocityWall * wallDistance + velocityFluid * (1 - wallDistance)
                  // Hence velocityWall = velocityFluid * (1 - 1/wallDistance)
                  util::Vector3D<double> velocityWall = util::Vector3D<double>(hydroVars.v_x,
                                                                               hydroVars.v_y,
                                                                               hydroVars.v_z)
                      * (1. - 1. / wallDistance);

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
                    // TODO I think we'll fail here if the next site out resides on a different core.

                    // Need some info about the next node away from the wall in this direction...
                    site_t nextSiteOutId = site.GetStreamedIndex(unstreamedDirection)
                        / D3Q15::NUMVECTORS;

                    if (log::Logger::ShouldDisplay<hemelb::log::Debug>())
                    {
                      if (nextSiteOutId < 0 || nextSiteOutId >= latDat->GetLocalFluidSiteCount())
                      {
                        log::Logger::Log<log::Debug, log::OnePerCore>("GZS "
                                                                      "boundary condition can't yet handle when the second fluid site away from "
                                                                      "a wall resides on a different core. The wall site was number %i on this core, "
                                                                      "the absent fluid site was in direction %i",
                                                                      siteIndex,
                                                                      unstreamedDirection);
                      }
                    }

                    geometry::Site nextSiteOut = latDat->GetSite(nextSiteOutId);

                    // Next, calculate its density, velocity and eqm distribution.
                    distribn_t nextNodeDensity, nextNodeV[3], nextNodeFEq[D3Q15::NUMVECTORS];
                    D3Q15::CalculateDensityVelocityFEq(nextSiteOut.GetFOld(),
                                                       nextNodeDensity,
                                                       nextNodeV[0],
                                                       nextNodeV[1],
                                                       nextNodeV[2],
                                                       nextNodeFEq);

                    // Obtain a second estimate, this time ignoring the fluid site closest to
                    // the wall. Interpolating the next site away and the site within the wall
                    // to the point on the wall itself (velocity 0):
                    // 0 = velocityWall * (1 + wallDistance) / 2 + velocityNextFluid * (1 - wallDistance)/2
                    // Rearranging gives velocityWall = velocityNextFluid * (wallDistance - 1)/(wallDistance+1)
                    util::Vector3D<double> velocityWallSecondEstimate =
                        util::Vector3D<double>(nextNodeV[0], nextNodeV[1], nextNodeV[2])
                            * (wallDistance - 1) / (wallDistance + 1);

                    // Next, we interpolate between the first and second estimates to improve the estimate.
                    // Extrapolate to obtain the velocity at the wall site.
                    for (int dimension = 0; dimension < 3; dimension++)
                    {
                      velocityWall[dimension] = wallDistance * velocityWall[dimension]
                          + (1. - wallDistance) * velocityWallSecondEstimate[dimension];
                    }

                    // Interpolate in the same way to get f_neq.
                    fNeqInUnstreamedDirection = wallDistance * fNeqInUnstreamedDirection
                        + (1. - wallDistance)
                            * (nextSiteOut.GetFOld()[unstreamedDirection]
                                - nextNodeFEq[unstreamedDirection]);
                  }

                  // Use a helper function to calculate the actual value of f_eq in the desired direction at the wall node.
                  // Note that we assume that the density is the same as at this node
                  distribn_t fEqTemp[D3Q15::NUMVECTORS];
                  D3Q15::CalculateFeq(hydroVars.density,
                                      velocityWall[0],
                                      velocityWall[1],
                                      velocityWall[2],
                                      fEqTemp);

                  // Collide and stream!
                  // TODO: It's not clear whether we should defer to the template collision type here
                  // or do a standard LBGK (implemented).
                  * (latDat->GetFNew(siteIndex * D3Q15::NUMVECTORS + unstreamedDirection)) =
                      fEqTemp[unstreamedDirection]
                          + (1.0 + lbmParams->GetOmega()) * fNeqInUnstreamedDirection;
                }
              }

              BaseStreamer<GuoZhengShi>::template UpdateMinsAndMaxes<tDoRayTracing>(hydroVars.v_x,
                                                                                    hydroVars.v_y,
                                                                                    hydroVars.v_z,
                                                                                    site,
                                                                                    hydroVars.GetFNeq().f,
                                                                                    hydroVars.density,
                                                                                    lbmParams,
                                                                                    visControl);
            }
          }

          template<bool tDoRayTracing>
          inline void DoPostStep(const site_t iFirstIndex,
                                 const site_t iSiteCount,
                                 const LbmParameters* iLbmParams,
                                 geometry::LatticeData* bLatDat,
                                 hemelb::vis::Control *iControl)
          {

          }

      };
    }
  }
}

#endif /* HEMELB_LB_STREAMERS_GUOZHENGSHI_H */
