#ifndef HEMELB_LB_COLLISIONS_IMPLEMENTATIONS_GUOZHENGSHI_H
#define HEMELB_LB_COLLISIONS_IMPLEMENTATIONS_GUOZHENGSHI_H

#include "lb/collisions/implementations/Implementation.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {
      namespace implementations
      {

        template<typename tCollisionOperator>
        class GuoZhengShi : public Implementation
        {

          public:
            template<bool tDoRayTracing>
            static void DoStreamAndCollide(WallCollision* mWallCollision,
                                           const site_t iFirstIndex,
                                           const site_t iSiteCount,
                                           const LbmParameters* iLbmParams,
                                           geometry::LatticeData* bLatDat,
                                           hemelb::vis::Control *iControl);

            template<bool tDoRayTracing>
            static void DoPostStep(WallCollision* mWallCollision,
                                   const site_t iFirstIndex,
                                   const site_t iSiteCount,
                                   const LbmParameters* iLbmParams,
                                   geometry::LatticeData* bLatDat,
                                   hemelb::vis::Control *iControl);

        };

        template<typename tCollisionOperator>
        template<bool tDoRayTracing>
        void GuoZhengShi<tCollisionOperator>::DoStreamAndCollide(WallCollision* mWallCollision,
                                                                 const site_t iFirstIndex,
                                                                 const site_t iSiteCount,
                                                                 const LbmParameters* iLbmParams,
                                                                 geometry::LatticeData* bLatDat,
                                                                 hemelb::vis::Control *iControl)
        {
          for (site_t lIndex = iFirstIndex; lIndex < (iFirstIndex + iSiteCount); lIndex++)
          {

            // First do a normal collision & streaming step, as if we were mid-fluid.
            // NOTE that we use the version that preserves f_old.
            // NOTE that this handily works out the equilibrium density, v_x, v_y and v_z for us
            distribn_t lFEq[15];
            distribn_t f_neq[15];
            distribn_t density, v_x, v_y, v_z;
            distribn_t* f = bLatDat->GetFOld(lIndex * D3Q15::NUMVECTORS);

            D3Q15::CalculateDensityVelocityFEq(f, density, v_x, v_y, v_z, lFEq);

            for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
            {
              * (bLatDat->GetFNew(bLatDat->GetStreamedIndex(lIndex, ii))) = f[ii]
                  + iLbmParams->Omega * (f_neq[ii] = f[ii] - lFEq[ii]);
            }

            // Now fill in the un-streamed-to distributions (those that point away from boundaries).
            for (unsigned int l = 1; l < D3Q15::NUMVECTORS; l++)
            {
              int lAwayFromWallIndex = D3Q15::INVERSEDIRECTIONS[l];

              if (bLatDat->HasBoundary(lIndex, l))
              {
                double delta = bLatDat->GetCutDistance(lIndex, l);
                double uWall[3];
                distribn_t fNeq;

                // Work out uw1 (noting that ub is 0 until we implement moving walls)
                uWall[0] = (1 - 1. / delta) * v_x;
                uWall[1] = (1 - 1. / delta) * v_y;
                uWall[2] = (1 - 1. / delta) * v_z;
                fNeq = f_neq[lAwayFromWallIndex];

                // Interpolate with uw2 if delta < 0.75
                if (delta < 0.75)
                {
                  // Only do the extra interpolation if there's gonna be a point there to interpolate from, i.e. there's no boundary
                  // in the direction of awayFromWallIndex
                  if (!bLatDat->HasBoundary(lIndex, lAwayFromWallIndex))
                  {
                    // Need some info about the next node away from the wall in this direction...
                    site_t nextIOut = bLatDat->GetStreamedIndex(lIndex, lAwayFromWallIndex)
                        / D3Q15::NUMVECTORS;
                    distribn_t nextNodeDensity, nextNodeV[3], nextNodeFEq[D3Q15::NUMVECTORS];

                    D3Q15::CalculateDensityVelocityFEq(bLatDat->GetFOld(nextIOut
                                                           * D3Q15::NUMVECTORS),
                                                       nextNodeDensity,
                                                       nextNodeV[0],
                                                       nextNodeV[1],
                                                       nextNodeV[2],
                                                       nextNodeFEq);

                    for (int a = 0; a < 3; a++)
                      uWall[a] = delta * uWall[a] + (1. - delta) * (delta - 1.) * nextNodeV[a]
                          / (1. + delta);

                    fNeq = delta * fNeq + (1. - delta)
                        * (*bLatDat->GetFOld(nextIOut * D3Q15::NUMVECTORS + lAwayFromWallIndex)
                            - nextNodeFEq[lAwayFromWallIndex]);
                  }
                  // If there's nothing to extrapolate from we, very lamely, do a 0VE-style operation to fill in the missing velocity.
                  else
                  {
                    for (int a = 0; a < 3; a++)
                      uWall[a] = 0.0;//delta * uWall[a];

                    fNeq = 0.0;//delta * fNeq;
                  }
                }

                // Use a helper function to calculate the actual value of f_eq in the desired direction at the wall node.
                // Note that we assume that the density is the same as at this node
                distribn_t fEqTemp[D3Q15::NUMVECTORS];
                D3Q15::CalculateFeq(density, uWall[0], uWall[1], uWall[2], fEqTemp);

                // Collide and stream!
                * (bLatDat->GetFNew(lIndex * D3Q15::NUMVECTORS + lAwayFromWallIndex))
                    = fEqTemp[lAwayFromWallIndex] + (1.0 + iLbmParams->Omega) * fNeq;
              }
            }

            UpdateMinsAndMaxes<tDoRayTracing> (v_x,
                                               v_y,
                                               v_z,
                                               lIndex,
                                               f_neq,
                                               density,
                                               bLatDat,
                                               iLbmParams,
                                               iControl);
          }
        }

        template<typename tCollisionOperator>
        template<bool tDoRayTracing>
        void GuoZhengShi<tCollisionOperator>::DoPostStep(WallCollision* mWallCollision,
                                                         const site_t iFirstIndex,
                                                         const site_t iSiteCount,
                                                         const LbmParameters* iLbmParams,
                                                         geometry::LatticeData* bLatDat,
                                                         hemelb::vis::Control *iControl)
        {

        }

      }
    }
  }
}

#endif /* HEMELB_LB_COLLISIONS_IMPLEMENTATIONS_GUOZHENGSHI_H */
