#ifndef HEMELB_LB_STREAMERS_GUOZHENGSHI_H
#define HEMELB_LB_STREAMERS_GUOZHENGSHI_H

#include "lb/streamers/BaseStreamer.h"

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {
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
          inline void DoStreamAndCollide(const site_t iFirstIndex,
                                         const site_t iSiteCount,
                                         const LbmParameters* iLbmParams,
                                         geometry::LatticeData* bLatDat,
                                         hemelb::vis::Control *iControl)
          {
            for (site_t lIndex = iFirstIndex; lIndex < (iFirstIndex + iSiteCount); lIndex++)
            {
              // First do a normal collision & streaming step, as if we were mid-fluid.
              distribn_t* f = bLatDat->GetFOld(lIndex * D3Q15::NUMVECTORS);
              kernels::HydroVars<typename CollisionType::CKernel> hydroVars(f);

              collider.CalculatePreCollision(hydroVars, lIndex - iFirstIndex);

              collider.Collide(iLbmParams, hydroVars);

              for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
              {
                * (bLatDat->GetFNew(bLatDat->GetStreamedIndex(lIndex, ii))) =
                    hydroVars.GetFPostCollision()[ii];

              }

              // Now fill in the un-streamed-to distributions (those that point away from boundaries).
              for (unsigned int l = 1; l < D3Q15::NUMVECTORS; l++)
              {
                if (bLatDat->HasBoundary(lIndex, l))
                {
                  int lAwayFromWallIndex = D3Q15::INVERSEDIRECTIONS[l];

                  double delta = bLatDat->GetCutDistance(lIndex, l);

                  // Work out uw1 (noting that ub is 0 until we implement moving walls)
                  double uWall[3];
                  uWall[0] = (1 - 1. / delta) * hydroVars.v_x;
                  uWall[1] = (1 - 1. / delta) * hydroVars.v_y;
                  uWall[2] = (1 - 1. / delta) * hydroVars.v_z;
                  distribn_t fNeqAwayFromWall = hydroVars.GetFNeq().f[lAwayFromWallIndex];

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
                      {
                        uWall[a] = delta * uWall[a]
                            - (1. - delta) * (1. - delta) * nextNodeV[a] / (1. + delta);
                      }

                      fNeqAwayFromWall = delta * fNeqAwayFromWall
                          + (1. - delta)
                              * (*bLatDat->GetFOld(nextIOut * D3Q15::NUMVECTORS
                                  + lAwayFromWallIndex) - nextNodeFEq[lAwayFromWallIndex]);
                    }
                    // If there's nothing to extrapolate from we, very lamely, do a 0VE-style operation to fill in the missing velocity.
                    else
                    {
                      for (int a = 0; a < 3; a++)
                      {
                        uWall[a] = 0.0; //delta * uWall[a];
                      }

                      fNeqAwayFromWall = 0.0; //delta * fNeq;
                    }
                  }

                  // Use a helper function to calculate the actual value of f_eq in the desired direction at the wall node.
                  // Note that we assume that the density is the same as at this node
                  distribn_t fEqTemp[D3Q15::NUMVECTORS];
                  D3Q15::CalculateFeq(hydroVars.density, uWall[0], uWall[1], uWall[2], fEqTemp);

                  // Collide and stream!
                  * (bLatDat->GetFNew(lIndex * D3Q15::NUMVECTORS + lAwayFromWallIndex)) =
                      fEqTemp[lAwayFromWallIndex]
                          + (1.0 + iLbmParams->GetOmega()) * fNeqAwayFromWall;
                }
              }

              BaseStreamer<GuoZhengShi>::template UpdateMinsAndMaxes<tDoRayTracing>(hydroVars.v_x,
                                                                                    hydroVars.v_y,
                                                                                    hydroVars.v_z,
                                                                                    lIndex,
                                                                                    hydroVars.GetFNeq().f,
                                                                                    hydroVars.density,
                                                                                    bLatDat,
                                                                                    iLbmParams,
                                                                                    iControl);
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
