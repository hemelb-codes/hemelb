#ifndef HEMELB_LB_COLLISIONS_STREAMANDCOLLIDE_H
#define HEMELB_LB_COLLISIONS_STREAMANDCOLLIDE_H

#include "lb/collisions/CollisionVisitor.h"
#include "lb/collisions/HFunction.h"
#include "util/utilityFunctions.h"
#include "D3Q15.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      template<bool tDoEntropic>
      class StreamAndCollide : public CollisionVisitor
      {
        public:
          virtual void VisitInletOutlet(InletOutletCollision* mInletOutletCollision,
                                        const bool iDoRayTracing,
                                        const site_t iFirstIndex,
                                        const site_t iSiteCount,
                                        const LbmParameters* iLbmParams,
                                        geometry::LatticeData* bLatDat,
                                        hemelb::vis::Control *iControl);

          virtual void VisitInletOutletWall(InletOutletWallCollision* mInletOutletWallCollision,
                                            const bool iDoRayTracing,
                                            const site_t iFirstIndex,
                                            const site_t iSiteCount,
                                            const LbmParameters* iLbmParams,
                                            geometry::LatticeData* bLatDat,
                                            hemelb::vis::Control *iControl);

          virtual void VisitMidFluid(MidFluidCollision* mMidFluidCollision,
                                     const bool iDoRayTracing,
                                     const site_t iFirstIndex,
                                     const site_t iSiteCount,
                                     const LbmParameters* iLbmParams,
                                     geometry::LatticeData* bLatDat,
                                     hemelb::vis::Control *iControl);

          virtual void VisitWall(WallCollision* mWallCollision,
                                 const bool iDoRayTracing,
                                 const site_t iFirstIndex,
                                 const site_t iSiteCount,
                                 const LbmParameters* iLbmParams,
                                 geometry::LatticeData* bLatDat,
                                 hemelb::vis::Control *iControl);

        private:
          template<bool tDoRayTracing>
          static void FInterpolation(WallCollision* mWallCollision,
                                     const site_t iFirstIndex,
                                     const site_t iSiteCount,
                                     const LbmParameters* iLbmParams,
                                     geometry::LatticeData* bLatDat,
                                     hemelb::vis::Control *iControl);

          template<bool tDoRayTracing>
          static void GuoZhengShi(WallCollision* mWallCollision,
                                  const site_t iFirstIndex,
                                  const site_t iSiteCount,
                                  const LbmParameters* iLbmParams,
                                  geometry::LatticeData* bLatDat,
                                  hemelb::vis::Control *iControl);

          template<bool tDoRayTracing>
          static void NonZeroVelocityBoundaryDensity(InletOutletCollision* mInletOutletCollision,
                                                     const site_t iFirstIndex,
                                                     const site_t iSiteCount,
                                                     const LbmParameters* iLbmParams,
                                                     geometry::LatticeData* bLatDat,
                                                     hemelb::vis::Control *iControl);

          template<bool tDoRayTracing>
          static void Regularised(WallCollision* mWallCollision,
                                  const site_t iFirstIndex,
                                  const site_t iSiteCount,
                                  const LbmParameters* iLbmParams,
                                  geometry::LatticeData* bLatDat,
                                  hemelb::vis::Control *iControl);

          template<bool tDoRayTracing>
          static void SimpleBounceBack(WallCollision* mWallCollision,
                                       const site_t iFirstIndex,
                                       const site_t iSiteCount,
                                       const LbmParameters* iLbmParams,
                                       geometry::LatticeData* bLatDat,
                                       hemelb::vis::Control *iControl);

          template<bool tDoRayTracing>
          static void SimpleCollideAndStream(MidFluidCollision* mMidFluidCollision,
                                             const site_t iFirstIndex,
                                             const site_t iSiteCount,
                                             const LbmParameters* iLbmParams,
                                             geometry::LatticeData* bLatDat,
                                             hemelb::vis::Control *iControl);

          template<bool tDoRayTracing>
          static void
          ZeroVelocityBoundaryDensity(InletOutletWallCollision* mInletOutletWallCollision,
                                      const site_t iFirstIndex,
                                      const site_t iSiteCount,
                                      const LbmParameters* iLbmParams,
                                      geometry::LatticeData* bLatDat,
                                      hemelb::vis::Control *iControl);

          template<bool tDoRayTracing>
          static void ZeroVelocityEquilibrium(WallCollision* mWallCollision,
                                              const site_t iFirstIndex,
                                              const site_t iSiteCount,
                                              const LbmParameters* iLbmParams,
                                              geometry::LatticeData* bLatDat,
                                              hemelb::vis::Control *iControl);

          static double getAlpha(distribn_t* lFOld, distribn_t* lFEq);

      }; /* End of StreamAndCollide definition */

      template<bool tDoEntropic>
      void StreamAndCollide<tDoEntropic>::VisitInletOutlet(InletOutletCollision* mInletOutletCollision,
                                                           const bool iDoRayTracing,
                                                           const site_t iFirstIndex,
                                                           const site_t iSiteCount,
                                                           const LbmParameters* iLbmParams,
                                                           geometry::LatticeData* bLatDat,
                                                           hemelb::vis::Control *iControl)
      {
        if (iDoRayTracing)
        {
          NonZeroVelocityBoundaryDensity<true> (mInletOutletCollision,
                                                iFirstIndex,
                                                iSiteCount,
                                                iLbmParams,
                                                bLatDat,
                                                iControl);
        }
        else
        {
          NonZeroVelocityBoundaryDensity<false> (mInletOutletCollision,
                                                 iFirstIndex,
                                                 iSiteCount,
                                                 iLbmParams,
                                                 bLatDat,
                                                 iControl);
        }
      }

      template<bool tDoEntropic>
      void StreamAndCollide<tDoEntropic>::VisitInletOutletWall(InletOutletWallCollision* mInletOutletWallCollision,
                                                               const bool iDoRayTracing,
                                                               const site_t iFirstIndex,
                                                               const site_t iSiteCount,
                                                               const LbmParameters* iLbmParams,
                                                               geometry::LatticeData* bLatDat,
                                                               hemelb::vis::Control *iControl)
      {
        if (iDoRayTracing)
        {
          ZeroVelocityBoundaryDensity<true> (mInletOutletWallCollision,
                                             iFirstIndex,
                                             iSiteCount,
                                             iLbmParams,
                                             bLatDat,
                                             iControl);
        }
        else
        {
          ZeroVelocityBoundaryDensity<false> (mInletOutletWallCollision,
                                              iFirstIndex,
                                              iSiteCount,
                                              iLbmParams,
                                              bLatDat,
                                              iControl);
        }
      }

      template<bool tDoEntropic>
      void StreamAndCollide<tDoEntropic>::VisitMidFluid(MidFluidCollision* mMidFluidCollision,
                                                        const bool iDoRayTracing,
                                                        const site_t iFirstIndex,
                                                        const site_t iSiteCount,
                                                        const LbmParameters* iLbmParams,
                                                        geometry::LatticeData* bLatDat,
                                                        hemelb::vis::Control *iControl)
      {
        if (iDoRayTracing)
        {
          SimpleCollideAndStream<true> (mMidFluidCollision,
                                        iFirstIndex,
                                        iSiteCount,
                                        iLbmParams,
                                        bLatDat,
                                        iControl);
        }
        else
        {
          SimpleCollideAndStream<false> (mMidFluidCollision,
                                         iFirstIndex,
                                         iSiteCount,
                                         iLbmParams,
                                         bLatDat,
                                         iControl);
        }
      }

      template<bool tDoEntropic>
      void StreamAndCollide<tDoEntropic>::VisitWall(WallCollision* mWallCollision,
                                                    const bool iDoRayTracing,
                                                    const site_t iFirstIndex,
                                                    const site_t iSiteCount,
                                                    const LbmParameters* iLbmParams,
                                                    geometry::LatticeData* bLatDat,
                                                    hemelb::vis::Control *iControl)
      {
        if (iDoRayTracing)
        {
          ZeroVelocityEquilibrium<true> (mWallCollision,
                                         iFirstIndex,
                                         iSiteCount,
                                         iLbmParams,
                                         bLatDat,
                                         iControl);
        }
        else
        {
          ZeroVelocityEquilibrium<false> (mWallCollision,
                                          iFirstIndex,
                                          iSiteCount,
                                          iLbmParams,
                                          bLatDat,
                                          iControl);
        }
      }

      template<bool tDoEntropic>
      template<bool tDoRayTracing>
      void StreamAndCollide<tDoEntropic>::FInterpolation(WallCollision* mWallCollision,
                                                         const site_t iFirstIndex,
                                                         const site_t iSiteCount,
                                                         const LbmParameters* iLbmParams,
                                                         geometry::LatticeData* bLatDat,
                                                         hemelb::vis::Control *iControl)
      {
        for (site_t lIndex = iFirstIndex; lIndex < (iFirstIndex + iSiteCount); lIndex++)
        {
          distribn_t* f = bLatDat->GetFOld(lIndex * D3Q15::NUMVECTORS);
          distribn_t density, v_x, v_y, v_z, f_neq[15];
          // Temporarily store f_eq in f_neq. Rectified later.
          D3Q15::CalculateDensityVelocityFEq(f, density, v_x, v_y, v_z, f_neq);

          for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
          {
            * (bLatDat->GetFNew(bLatDat->GetStreamedIndex(lIndex, ii))) = f[ii]
                += iLbmParams->Omega * (f_neq[ii] = f[ii] - f_neq[ii]);
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

      template<bool tDoEntropic>
      template<bool tDoRayTracing>
      void StreamAndCollide<tDoEntropic>::GuoZhengShi(WallCollision* mWallCollision,
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
            * (bLatDat->GetFNew(bLatDat->GetStreamedIndex(lIndex, ii))) = f[ii] + iLbmParams->Omega
                * (f_neq[ii] = f[ii] - lFEq[ii]);
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

                  D3Q15::CalculateDensityVelocityFEq(bLatDat->GetFOld(nextIOut * D3Q15::NUMVECTORS),
                                                     nextNodeDensity,
                                                     nextNodeV[0],
                                                     nextNodeV[1],
                                                     nextNodeV[2],
                                                     nextNodeFEq);

                  for (int a = 0; a < 3; a++)
                    uWall[a] = delta * uWall[a] + (1. - delta) * (delta - 1.) * nextNodeV[a] / (1.
                        + delta);

                  fNeq = delta * fNeq + (1. - delta) * (*bLatDat->GetFOld(nextIOut
                      * D3Q15::NUMVECTORS + lAwayFromWallIndex) - nextNodeFEq[lAwayFromWallIndex]);
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

      template<bool tDoEntropic>
      template<bool tDoRayTracing>
      void StreamAndCollide<tDoEntropic>::NonZeroVelocityBoundaryDensity(InletOutletCollision* mInletOutletCollision,
                                                                         const site_t iFirstIndex,
                                                                         const site_t iSiteCount,
                                                                         const LbmParameters* iLbmParams,
                                                                         geometry::LatticeData* bLatDat,
                                                                         hemelb::vis::Control *iControl)
      {
        for (site_t lIndex = iFirstIndex; lIndex < (iFirstIndex + iSiteCount); lIndex++)
        {
          distribn_t* lFOld = bLatDat->GetFOld(lIndex * D3Q15::NUMVECTORS);
          distribn_t lFNeq[15];
          distribn_t lVx, lVy, lVz, lDummyDensity, lDensity;

          lDensity
              = (*mInletOutletCollision).getBoundaryDensityArray(bLatDat->GetBoundaryId(lIndex));

          D3Q15::CalculateDensityAndVelocity(lFOld, lDummyDensity, lVx, lVy, lVz);

          for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
          {
            lFNeq[ii] = lFOld[ii];
          }

          // Temporarily store FEq in lFNeq (rectified later).
          if (tDoEntropic)
            D3Q15::CalculateEntropicFeq(lDensity, lVx, lVy, lVz, lFOld);
          else
            D3Q15::CalculateFeq(lDensity, lVx, lVy, lVz, lFOld);

          for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
          {
            * (bLatDat->GetFNew(bLatDat->GetStreamedIndex(lIndex, ii))) = lFOld[ii];
          }

          for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
          {
            lFNeq[ii] -= lFOld[ii];
          }

          UpdateMinsAndMaxes<tDoRayTracing> (lVx,
                                             lVy,
                                             lVz,
                                             lIndex,
                                             lFNeq,
                                             lDensity,
                                             bLatDat,
                                             iLbmParams,
                                             iControl);
        }
      }

      template<bool tDoEntropic>
      template<bool tDoRayTracing>
      void StreamAndCollide<tDoEntropic>::Regularised(WallCollision* mWallCollision,
                                                      const site_t iFirstIndex,
                                                      const site_t iSiteCount,
                                                      const LbmParameters* iLbmParams,
                                                      geometry::LatticeData* bLatDat,
                                                      hemelb::vis::Control *iControl)
      {
        for (site_t lIndex = iFirstIndex; lIndex < (iFirstIndex + iSiteCount); lIndex++)
        {
          distribn_t *f = bLatDat->GetFOld(lIndex * D3Q15::NUMVECTORS);
          distribn_t density, v_x, v_y, v_z;
          distribn_t f_neq[15];

          // First calculate the density and macro-velocity
          // TEMPORARILY STORE f_eq IN f_neq BUT THE FUNCTION RETURNS f_eq. THIS IS SORTED
          // OUT IN A SUBSEQUENT FOR LOOP.
          D3Q15::CalculateDensityVelocityFEq(f, density, v_x, v_y, v_z, f_neq);

          // To evaluate PI, first let unknown particle populations take value given by bounce-back of off-equilibrium parts
          // (fi = fiEq + fopp(i) - fopp(i)Eq)
          distribn_t fTemp[15];

          for (int l = 0; l < 15; l++)
            fTemp[l] = f[l];

          if (bLatDat->HasBoundary(lIndex, 1))
          {
            fTemp[1] = f[2] + (2.0 / 3.0) * v_x;
          }
          if (bLatDat->HasBoundary(lIndex, 2))
          {
            fTemp[2] = f[1] - (2.0 / 3.0) * v_x;
          }
          if (bLatDat->HasBoundary(lIndex, 3))
          {
            fTemp[3] = f[4] + (2.0 / 3.0) * v_y;
          }
          if (bLatDat->HasBoundary(lIndex, 4))
          {
            fTemp[4] = f[3] - (2.0 / 3.0) * v_y;
          }
          if (bLatDat->HasBoundary(lIndex, 5))
          {
            fTemp[5] = f[6] + (2.0 / 3.0) * v_z;
          }
          if (bLatDat->HasBoundary(lIndex, 6))
          {
            fTemp[6] = f[5] - (2.0 / 3.0) * v_z;
          }
          if (bLatDat->HasBoundary(lIndex, 7))
          {
            fTemp[7] = f[8] + (2.0 / 24.0) * ( (v_x + v_y) + v_z);
          }
          if (bLatDat->HasBoundary(lIndex, 8))
          {
            fTemp[8] = f[7] - (2.0 / 24.0) * ( (v_x + v_y) + v_z);
          }
          if (bLatDat->HasBoundary(lIndex, 9))
          {
            fTemp[9] = f[10] + (2.0 / 24.0) * ( (v_x + v_y) - v_z);
          }
          if (bLatDat->HasBoundary(lIndex, 10))
          {
            fTemp[10] = f[9] - (2.0 / 24.0) * ( (v_x + v_y) - v_z);
          }
          if (bLatDat->HasBoundary(lIndex, 11))
          {
            fTemp[11] = f[12] + (2.0 / 24.0) * ( (v_x - v_y) + v_z);
          }
          if (bLatDat->HasBoundary(lIndex, 12))
          {
            fTemp[12] = f[11] - (2.0 / 24.0) * ( (v_x - v_y) + v_z);
          }
          if (bLatDat->HasBoundary(lIndex, 13))
          {
            fTemp[13] = f[14] + (2.0 / 24.0) * ( (v_x - v_y) - v_z);
          }
          if (bLatDat->HasBoundary(lIndex, 14))
          {
            fTemp[14] = f[13] - (2.0 / 24.0) * ( (v_x - v_y) - v_z);
          }

          // UP TO THIS POINT, F_NEQ ACTUALLY CONTAINS F_EQ. AT THIS
          // STAGE WE REPLACE IT WITH THE ACTUAL NON-EQ VALUE, POST
          // BOUNCING_BACK WHERE NEEDED.
          for (int l = 0; l < 15; l++)
          {
            f_neq[l] = fTemp[l] - f_neq[l];
          }

          distribn_t density_1 = 1. / density;
          distribn_t v_xx = v_x * v_x;
          distribn_t v_yy = v_y * v_y;
          distribn_t v_zz = v_z * v_z;

          // PI = sum_i e_i e_i f_i
          distribn_t piMatrix[3][3];

          distribn_t diagSum = f_neq[7] + f_neq[8] + f_neq[9] + f_neq[10] + f_neq[11] + f_neq[12]
              + f_neq[13] + f_neq[14];

          piMatrix[0][0] = f_neq[1] + f_neq[2] + diagSum;
          piMatrix[0][1] = f_neq[7] + f_neq[8] + f_neq[9] + f_neq[10] - (f_neq[11] + f_neq[12]
              + f_neq[13] + f_neq[14]);
          piMatrix[0][2] = f_neq[7] + f_neq[8] + f_neq[11] + f_neq[12] - (f_neq[9] + f_neq[10]
              + f_neq[13] + f_neq[14]);
          piMatrix[1][0] = piMatrix[0][1];
          piMatrix[1][1] = f_neq[3] + f_neq[4] + diagSum;
          piMatrix[1][2] = f_neq[7] + f_neq[8] + f_neq[13] + f_neq[14] - (f_neq[9] + f_neq[10]
              + f_neq[11] + f_neq[12]);
          piMatrix[2][0] = piMatrix[0][2];
          piMatrix[2][1] = piMatrix[1][2];
          piMatrix[2][2] = f_neq[5] + f_neq[6] + diagSum;

          for (int m = 0; m < 3; m++)
            for (int n = 0; n < 3; n++)
              piMatrix[m][n] /= (2.0 * Cs2 * Cs2);

          // Qi = e_i e_i - (speed of sound ^ 2) * Identity
          // Then gi = fiEq + t_i (the 2/9, 1/9, 1/72 stuff) (Qi . PI (inner product)) / 2 * speed of sound^4
          // Or:  gi = fiEq + t_i (the 2/9, 1/9, 1/72 stuff) ((e_i e_i . PI (inner product)) / 2 * speed of sound^4 - specialNumber)
          distribn_t specialNumber = (2.0 / 9.0) * Cs2 * (piMatrix[0][0] + piMatrix[1][1]
              + piMatrix[2][2]);
          distribn_t piMatrixSum = piMatrix[0][0] + piMatrix[0][1] + piMatrix[0][2]
              + piMatrix[1][0] + piMatrix[1][1] + piMatrix[1][2] + piMatrix[2][0] + piMatrix[2][1]
              + piMatrix[2][2];

          // The gi (here; f) are then collided and streamed
          * (bLatDat->GetFNew(bLatDat->GetStreamedIndex(lIndex, 0))) = ( (2.0 / 9.0) * density
              - (1.0 / 3.0) * ( (v_xx + v_yy + v_zz) * density_1)) + (1.0 + iLbmParams->Omega)
              * (f_neq[0] = -specialNumber);

          distribn_t temp1 = (1.0 / 9.0) * density - (1.0 / 6.0) * ( (v_xx + v_yy + v_zz)
              * density_1);
          specialNumber *= 1.0 / 2.0;

          // Now apply bounce-back to the components that require it, from fTemp
          site_t lStreamTo[15];
          for (int l = 1; l < 15; l++)
          {
            if (bLatDat->HasBoundary(lIndex, l))
            {
              lStreamTo[l] = lIndex * 15 + D3Q15::INVERSEDIRECTIONS[l];
            }
            else
            {
              lStreamTo[l] = bLatDat->GetStreamedIndex(lIndex, l);
            }
          }

          * (bLatDat->GetFNew(lStreamTo[1])) = temp1 + (0.5 * density_1) * v_xx + (1.0 / 3.0) * v_x
              + (1.0 + iLbmParams->Omega) * (f_neq[1] = (1.0 / 9.0) * piMatrix[0][0]
                  - specialNumber); // (+1, 0, 0)
          * (bLatDat->GetFNew(lStreamTo[2])) = temp1 + (0.5 * density_1) * v_xx - (1.0 / 3.0) * v_x
              + (1.0 + iLbmParams->Omega) * (f_neq[2] = (1.0 / 9.0) * piMatrix[0][0]
                  - specialNumber); // (+1, 0, 0)

          * (bLatDat->GetFNew(lStreamTo[3])) = temp1 + (0.5 * density_1) * v_yy + (1.0 / 3.0) * v_y
              + (1.0 + iLbmParams->Omega) * (f_neq[3] = (1.0 / 9.0) * piMatrix[1][1]
                  - specialNumber); // (0, +1, 0)
          * (bLatDat->GetFNew(lStreamTo[4])) = temp1 + (0.5 * density_1) * v_yy - (1.0 / 3.0) * v_y
              + (1.0 + iLbmParams->Omega) * (f_neq[4] = (1.0 / 9.0) * piMatrix[1][1]
                  - specialNumber); // (0, +1, 0)

          * (bLatDat->GetFNew(lStreamTo[5])) = temp1 + (0.5 * density_1) * v_zz + (1.0 / 3.0) * v_z
              + (1.0 + iLbmParams->Omega) * (f_neq[5] = (1.0 / 9.0) * piMatrix[2][2]
                  - specialNumber); // (0, +1, 0)
          * (bLatDat->GetFNew(lStreamTo[6])) = temp1 + (0.5 * density_1) * v_zz - (1.0 / 3.0) * v_z
              + (1.0 + iLbmParams->Omega) * (f_neq[6] = (1.0 / 9.0) * piMatrix[2][2]
                  - specialNumber); // (0, +1, 0)

          temp1 *= (1.0 / 8.0);
          specialNumber *= (1.0 / 8.0);

          distribn_t temp2 = (v_x + v_y) + v_z;

          * (bLatDat->GetFNew(lStreamTo[7])) = temp1 + (1.0 / 16.0) * density_1 * temp2 * temp2
              + (1.0 / 24.0) * temp2 + (1.0 + iLbmParams->Omega) * (f_neq[7] = ( (1.0 / 72.0)
              * piMatrixSum - specialNumber)); // (+1, +1, +1)
          * (bLatDat->GetFNew(lStreamTo[8])) = temp1 + (1.0 / 16.0) * density_1 * temp2 * temp2
              + (-1.0 / 24.0) * temp2 + (1.0 + iLbmParams->Omega) * (f_neq[8] = ( (1.0 / 72.0)
              * piMatrixSum - specialNumber)); // (-1, -1, -1)

          temp2 = (v_x + v_y) - v_z;

          * (bLatDat->GetFNew(lStreamTo[9])) = temp1 + (1.0 / 16.0) * density_1 * temp2 * temp2
              + (1.0 / 24.0) * temp2 + (1.0 + iLbmParams->Omega) * (f_neq[9] = ( (1.0 / 72.0)
              * (piMatrixSum - 4.0 * (piMatrix[0][2] + piMatrix[1][2])) - specialNumber)); // (+1, +1, -1)
          * (bLatDat->GetFNew(lStreamTo[10])) = temp1 + (1.0 / 16.0) * density_1 * temp2 * temp2
              + (-1.0 / 24.0) * temp2 + (1.0 + iLbmParams->Omega) * (f_neq[10] = ( (1.0 / 72.0)
              * (piMatrixSum - 4.0 * (piMatrix[0][2] + piMatrix[1][2])) - specialNumber)); // (-1, -1, +1)

          temp2 = (v_x - v_y) + v_z;

          * (bLatDat->GetFNew(lStreamTo[11])) = temp1 + (1.0 / 16.0) * density_1 * temp2 * temp2
              + (1.0 / 24.0) * temp2 + (1.0 + iLbmParams->Omega) * (f_neq[11] = ( (1.0 / 72.0)
              * (piMatrixSum - 4.0 * (piMatrix[0][1] + piMatrix[1][2])) - specialNumber)); // (+1, -1, +1)
          * (bLatDat->GetFNew(lStreamTo[12])) = temp1 + (1.0 / 16.0) * density_1 * temp2 * temp2
              + (-1.0 / 24.0) * temp2 + (1.0 + iLbmParams->Omega) * (f_neq[12] = ( (1.0 / 72.0)
              * (piMatrixSum - 4.0 * (piMatrix[0][1] + piMatrix[1][2])) - specialNumber)); // (-1, +1, -1)

          temp2 = (v_x - v_y) - v_z;

          * (bLatDat->GetFNew(lStreamTo[13])) = temp1 + (1.0 / 16.0) * density_1 * temp2 * temp2
              + (1.0 / 24.0) * temp2 + (1.0 + iLbmParams->Omega) * (f_neq[13] = ( (1.0 / 72.0)
              * (piMatrixSum - 4.0 * (piMatrix[0][1] + piMatrix[0][2])) - specialNumber)); // (+1, -1, -1)
          * (bLatDat->GetFNew(lStreamTo[14])) = temp1 + (1.0 / 16.0) * density_1 * temp2 * temp2
              + (-1.0 / 24.0) * temp2 + (1.0 + iLbmParams->Omega) * (f_neq[14] = ( (1.0 / 72.0)
              * (piMatrixSum - 4.0 * (piMatrix[0][1] + piMatrix[0][2])) - specialNumber)); // (-1, +1, +1)

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

      template<bool tDoEntropic>
      template<bool tDoRayTracing>
      void StreamAndCollide<tDoEntropic>::SimpleBounceBack(WallCollision* mWallCollision,
                                                           const site_t iFirstIndex,
                                                           const site_t iSiteCount,
                                                           const LbmParameters* iLbmParams,
                                                           geometry::LatticeData* bLatDat,
                                                           hemelb::vis::Control *iControl)
      {
        for (site_t lIndex = iFirstIndex; lIndex < (iFirstIndex + iSiteCount); lIndex++)
        {
          distribn_t *lFOld = bLatDat->GetFOld(lIndex * D3Q15::NUMVECTORS);
          distribn_t lFNeq[D3Q15::NUMVECTORS];
          distribn_t lVx, lVy, lVz, lDensity;

          for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
          {
            lFNeq[ii] = lFOld[ii];
          }

          D3Q15::CalculateDensityVelocityFEq(lFNeq, lDensity, lVx, lVy, lVz, lFOld);

          for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
          {
            // The actual bounce-back lines, including streaming and collision. Basically swap the non-equilibrium components of f in each of the opposing pairs of directions.
            site_t lStreamTo = (bLatDat->HasBoundary(lIndex, ii))
              ? ( (lIndex * D3Q15::NUMVECTORS) + D3Q15::INVERSEDIRECTIONS[ii])
              : bLatDat->GetStreamedIndex(lIndex, ii);

            // Remember, oFNeq currently hold the equilibrium distribution. We
            // simultaneously use this and correct it, here.
            * (bLatDat->GetFNew(lStreamTo)) = lFOld[ii] += iLbmParams->Omega * (lFNeq[ii]
                -= lFOld[ii]);
          }

          UpdateMinsAndMaxes<tDoRayTracing> (lVx,
                                             lVy,
                                             lVz,
                                             lIndex,
                                             lFNeq,
                                             lDensity,
                                             bLatDat,
                                             iLbmParams,
                                             iControl);
        }
      }

      template<bool tDoEntropic>
      template<bool tDoRayTracing>
      void StreamAndCollide<tDoEntropic>::SimpleCollideAndStream(MidFluidCollision* mMidFluidCollision,
                                                                 const site_t iFirstIndex,
                                                                 const site_t iSiteCount,
                                                                 const LbmParameters* iLbmParams,
                                                                 geometry::LatticeData* bLatDat,
                                                                 hemelb::vis::Control *iControl)
      {
        for (site_t iIndex = iFirstIndex; iIndex < (iFirstIndex + iSiteCount); iIndex++)
        {
          distribn_t* lFOld = bLatDat->GetFOld(iIndex * D3Q15::NUMVECTORS);
          distribn_t lDensity, lVx, lVy, lVz;
          distribn_t lFNeq[D3Q15::NUMVECTORS];
          double alpha;

          // Temporarily store f_eq in f_neq (rectified in next statement)
          if (tDoEntropic)
          {
            D3Q15::CalculateEntropicDensityVelocityFEq(lFOld, lDensity, lVx, lVy, lVz, lFNeq);
            alpha = getAlpha(lFOld, lFNeq); // lFNeq is actually lFeq at the moment
          }
          else
          {
            D3Q15::CalculateDensityVelocityFEq(lFOld, lDensity, lVx, lVy, lVz, lFNeq);
          }

          for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
          {
            if (tDoEntropic)
            {
              * (bLatDat->GetFNew(bLatDat->GetStreamedIndex(iIndex, ii))) = lFOld[ii] += alpha
                  * iLbmParams->Beta * (lFNeq[ii] = lFOld[ii] - lFNeq[ii]);
            }
            else
            {
              * (bLatDat->GetFNew(bLatDat->GetStreamedIndex(iIndex, ii))) = lFOld[ii]
                  += iLbmParams->Omega * (lFNeq[ii] = lFOld[ii] - lFNeq[ii]);
            }
          }

          UpdateMinsAndMaxes<tDoRayTracing> (lVx,
                                             lVy,
                                             lVz,
                                             iIndex,
                                             lFNeq,
                                             lDensity,
                                             bLatDat,
                                             iLbmParams,
                                             iControl);
        }
      }

      template<bool tDoEntropic>
      template<bool tDoRayTracing>
      void StreamAndCollide<tDoEntropic>::ZeroVelocityBoundaryDensity(InletOutletWallCollision* mInletOutletWallCollision,
                                                                      const site_t iFirstIndex,
                                                                      const site_t iSiteCount,
                                                                      const LbmParameters* iLbmParams,
                                                                      geometry::LatticeData* bLatDat,
                                                                      hemelb::vis::Control *iControl)
      {
        for (site_t iIndex = iFirstIndex; iIndex < (iFirstIndex + iSiteCount); iIndex++)
        {
          distribn_t* lFOld = bLatDat->GetFOld(iIndex * D3Q15::NUMVECTORS);
          distribn_t lFNeq[D3Q15::NUMVECTORS];
          distribn_t lDensity;

          lDensity
              = (*mInletOutletWallCollision).getBoundaryDensityArray(bLatDat->GetBoundaryId(iIndex));

          for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
          {
            lFNeq[ii] = lFOld[ii];
          }

          // Temporarily store FEq in FNeq
          if (tDoEntropic)
            D3Q15::CalculateEntropicFeq(lDensity, 0.0, 0.0, 0.0, lFOld);
          else
            D3Q15::CalculateFeq(lDensity, 0.0, 0.0, 0.0, lFOld);

          for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
          {
            * (bLatDat->GetFNew(bLatDat->GetStreamedIndex(iIndex, ii))) = lFOld[ii];
            lFNeq[ii] -= lFOld[ii];
          }

          UpdateMinsAndMaxes<tDoRayTracing> (0.0,
                                             0.0,
                                             0.0,
                                             iIndex,
                                             lFNeq,
                                             lDensity,
                                             bLatDat,
                                             iLbmParams,
                                             iControl);
        }
      }

      template<bool tDoEntropic>
      template<bool tDoRayTracing>
      void StreamAndCollide<tDoEntropic>::ZeroVelocityEquilibrium(WallCollision* mWallCollision,
                                                                  const site_t iFirstIndex,
                                                                  const site_t iSiteCount,
                                                                  const LbmParameters* iLbmParams,
                                                                  geometry::LatticeData* bLatDat,
                                                                  hemelb::vis::Control *iControl)
      {
        for (site_t lIndex = iFirstIndex; lIndex < (iFirstIndex + iSiteCount); lIndex++)
        {
          distribn_t* lFOld = bLatDat->GetFOld(lIndex * D3Q15::NUMVECTORS);
          distribn_t lFNeq[D3Q15::NUMVECTORS];
          distribn_t lDensity;

          lDensity = 0.0;

          for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
          {
            lDensity += lFOld[ii];
          }

          for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
          {
            lFNeq[ii] = lFOld[ii];
          }

          // Temporarily store FEq in lFNeq
          if (tDoEntropic)
            D3Q15::CalculateEntropicFeq(lDensity, 0.0, 0.0, 0.0, lFOld);
          else
            D3Q15::CalculateFeq(lDensity, 0.0, 0.0, 0.0, lFOld);

          for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
          {
            * (bLatDat->GetFNew(bLatDat->GetStreamedIndex(lIndex, ii))) = lFOld[ii];
            lFNeq[ii] -= lFOld[ii];
          }

          UpdateMinsAndMaxes<tDoRayTracing> (0.0,
                                             0.0,
                                             0.0,
                                             lIndex,
                                             lFNeq,
                                             lDensity,
                                             bLatDat,
                                             iLbmParams,
                                             iControl);
        }
      }

      template<bool tDoEntropic>
      double StreamAndCollide<tDoEntropic>::getAlpha(distribn_t* lF, distribn_t* lFEq)
      {
        bool small = true;
        distribn_t deviation[D3Q15::NUMVECTORS];

        for (unsigned int i = 0; i < D3Q15::NUMVECTORS; i++)
        {
          if (fabs(/* deviation[i] = */ (lFEq[i] - lF[i]) / lF[i]) > 0.01)
          {
            small = false;
            break;
          }
        }

        if (small)
          return 2.0;

        /*
         if (small)
         {
         double a1 = 0.0, a2 = 0.0, a3 = 0.0, a4 = 0.0;
         double deviation2[D3Q15::NUMVECTORS]; //Will be holding deviation[i]^n * (f_eq - f)

         for (unsigned int i = 0; i < D3Q15::NUMVECTORS; i++)
         {
         a1 += (deviation2[i] = deviation[i] * (lFEq[i] - lF[i]));
         a2 += (deviation2[i] *= deviation[i]);
         a3 += (deviation2[i] *= deviation[i]);
         a4 += (deviation2[i] * deviation[i]);
         }

         a1 *= 0.5;
         a2 *= -1.0 / 6.0;
         a3 *= 1.0 / 12.0;
         a4 *= -0.05;

         // Alpha tends to 2 as deviation -> 0. The evaluation of alpha below fails if a1 is 0 so to prevent
         // NaNs and to save some computation we just return 2
         if (a1 < 1.0E-10)
         return 2.0;

         return (2.0 + 4.0 * (4 * a1 * a2 * (a2 + 5.0 * a3) - a1 * a1 * (a2 + 2.0 * a3
         + 4.0 * a4) - 20.0 * a2 * a2 * a2) / (a1 * a1 * a1));
         }
         */

        HFunction HFunc(lF, lFEq);

        return (hemelb::util::NumericalFunctions::NewtonRaphson(&HFunc, 2.0, 1.0E-6));

      }

    }
  }
}

#endif /* HEMELB_LB_COLLISIONS_STREAMANDCOLLIDE_H */
