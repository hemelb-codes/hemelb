
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_LBTESTS_VIRTUALSITEIOLETSTREAMERTESTS_H
#define HEMELB_UNITTESTS_LBTESTS_VIRTUALSITEIOLETSTREAMERTESTS_H

#include <cppunit/TestFixture.h>
#include <iostream>
#include <sstream>

#include "lb/streamers/Streamers.h"
#include "geometry/SiteData.h"
#include "util/utilityFunctions.h"
#include "unittests/helpers/FourCubeBasedTestFixture.h"

namespace hemelb
{
  namespace unittests
  {
    namespace lbtests
    {
      //static const distribn_t allowedError = 1e-10;
      /**
       * Helper class that exposes implementation details to the tests.
       */
      using lb::streamers::VSExtra;
      class LocalVSExtra : public VSExtra<lb::lattices::D3Q15>
      {
        public:
          //typedef hemelb::lb::streamers::VSExtra::UnitVec UnitVec;
          LocalVSExtra(lb::iolets::InOutLet& iolet) :
            VSExtra<lb::lattices::D3Q15> (iolet)
          {
          }
          const UnitVec& GetE1() const
          {
            return e1;
          }
          const UnitVec& GetE2() const
          {
            return e2;
          }
          const UnitVec& GetE3() const
          {
            return n;
          }
      };

      /**
       * VirtualSiteIoletStreamerTests:
       *
       * This class tests the streamer implementations. We assume the collision operators are
       * correct (as they're tested elsewhere), then compare the post-streamed values with
       * the values we expect to have been streamed there.
       */
      class VirtualSiteIoletStreamerTests : public helpers::FourCubeBasedTestFixture
      {
          CPPUNIT_TEST_SUITE( VirtualSiteIoletStreamerTests);
          CPPUNIT_TEST( TestVirtualSiteMatrixInverse);
          CPPUNIT_TEST( TestVirtualSiteConstruction);
          CPPUNIT_TEST( TestStreamerInitialisation);
          CPPUNIT_TEST( TestStep);CPPUNIT_TEST_SUITE_END();

        public:
          typedef lb::lattices::D3Q15 Lattice;
          typedef lb::kernels::LBGK<Lattice> Kernel;
          typedef lb::collisions::Normal<Kernel> Collision;
          typedef lb::streamers::VirtualSite<Lattice> VirtualSite;
          typedef lb::iolets::InOutLetCosine InOutLetCosine;
          typedef lb::streamers::RSHV RSHV;

          InOutLetCosine* GetIolet(lb::iolets::BoundaryValues* iolets)
          {
            InOutLetCosine* ans =
                dynamic_cast<lb::iolets::InOutLetCosine*> (iolets->GetLocalIolet(0));
            CPPUNIT_ASSERT(ans != NULL);
            return ans;
          }
          void setUp()
          {

            FourCubeBasedTestFixture::setUp();
            propertyCache = new lb::MacroscopicPropertyCache(*simState, *latDat);

            lattice = &Lattice::GetLatticeInfo();

            inletBoundary = new lb::iolets::BoundaryValues(geometry::INLET_TYPE,
                                                           latDat,
                                                           simConfig->GetInlets(),
                                                           simState,
                                                           Comms(),
                                                           *unitConverter);
            InOutLetCosine* inlet = GetIolet(inletBoundary);
            // We have to make the outlet sane and consistent with the geometry now.
            inlet->SetNormal(util::Vector3D<Dimensionless>(0, 0, 1));
            PhysicalPosition inletCentre(2.5, 2.5, 0.5);
            inletCentre *= simConfig->GetVoxelSize();
            inlet->SetPosition(unitConverter->ConvertPositionToLatticeUnits(inletCentre));
            // Want to set the density gradient to be 0.01 in lattice units,
            // starting at 1.0 at the outlet.
            inlet->SetPressureAmp(0.);
            inlet->SetPressureMean(unitConverter->ConvertPressureToPhysicalUnits(1.04 * Cs2));

            // Same for the cut distances of outlet sites.
            for (unsigned i = 1; i <= 4; ++i)
            {
              for (unsigned j = 1; j <= 4; ++j)
              {
                site_t siteId = latDat->GetContiguousSiteId(LatticeVector(i, j, 1));
                geometry::Site < geometry::LatticeData > site = latDat->GetSite(siteId);
                for (Direction p = 0; p < lattice->GetNumVectors(); ++p)
                {
                  if (lattice->GetVector(p).z < 0.)
                  {
                    // Sanity check
                    CPPUNIT_ASSERT(site.HasIolet(p));
                    // Set the cut distance to half way
                    latDat->SetBoundaryDistance(siteId, p, 0.5);
                  }
                }
              }
            }

            outletBoundary = new lb::iolets::BoundaryValues(geometry::OUTLET_TYPE,
                                                            latDat,
                                                            simConfig->GetOutlets(),
                                                            simState,
                                                            Comms(),
                                                            *unitConverter);

            InOutLetCosine* outlet = GetIolet(outletBoundary);
            // We have to make the outlet sane and consistent with the geometry now.
            outlet->SetNormal(util::Vector3D<Dimensionless>(0, 0, -1));
            PhysicalPosition outletCentre(2.5, 2.5, 4.5);
            outletCentre *= simConfig->GetVoxelSize();
            outlet->SetPosition(unitConverter->ConvertPositionToLatticeUnits(outletCentre));
            outlet->SetPressureAmp(0.);
            outlet->SetPressureMean(unitConverter->ConvertPressureToPhysicalUnits(1.0 * Cs2));

            // Same for the cut distances of outlet sites.
            for (unsigned i = 1; i <= 4; ++i)
            {
              for (unsigned j = 1; j <= 4; ++j)
              {
                site_t siteId = latDat->GetContiguousSiteId(LatticeVector(i, j, 4));
                geometry::Site < geometry::LatticeData > site = latDat->GetSite(siteId);
                for (Direction p = 0; p < lattice->GetNumVectors(); ++p)
                {
                  if (lattice->GetVector(p).z > 0.)
                  {
                    // Sanity check
                    CPPUNIT_ASSERT(site.HasIolet(p));
                    // Set the cut distance to half way
                    latDat->SetBoundaryDistance(siteId, p, 0.5);
                  }
                }
              }
            }
          }

          void tearDown()
          {
            delete propertyCache;
            delete inletBoundary;
            delete outletBoundary;
            FourCubeBasedTestFixture::tearDown();
          }

          void TestVirtualSiteMatrixInverse()
          {
            distribn_t inv[3][3];
            distribn_t det;

            // Try the identity matrix first
            distribn_t eye[3][3] = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
            det = VirtualSite::Matrix3DInverse(eye, inv);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(1., det, allowedError);
            for (unsigned i = 0; i < 3; ++i)
              for (unsigned j = 0; j < 3; ++j)
                CPPUNIT_ASSERT_DOUBLES_EQUAL(eye[i][j], inv[i][j], allowedError);

            // arbitrary matrix with nice ish values
            distribn_t a[3][3] = { { 1, 2, 3 }, { 0, 4, 5 }, { 1, 0, 6 } };
            distribn_t aInv22[3][3] = { { 24, -12, -2 }, { 5, 3, -5 }, { -4, 2, 4 } };

            det = VirtualSite::Matrix3DInverse(a, inv);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(22., det, allowedError);
            for (unsigned i = 0; i < 3; ++i)
              for (unsigned j = 0; j < 3; ++j)
                CPPUNIT_ASSERT_DOUBLES_EQUAL(aInv22[i][j] / 22, inv[i][j], allowedError);

          }

          void TestVirtualSiteConstruction()
          {
            InOutLetCosine& outlet = *GetIolet(outletBoundary);
            LocalVSExtra extra(outlet);
            for (unsigned i = 0; i < 3; ++i)
            {
              // Expected basis for iolet coords
              util::Vector3D<Dimensionless> e1(-1, 0, 0), e2(0, 1, 0), e3(0, 0, -1);
              CPPUNIT_ASSERT_DOUBLES_EQUAL(e1[i], extra.GetE1()[i], allowedError);
              CPPUNIT_ASSERT_DOUBLES_EQUAL(e2[i], extra.GetE2()[i], allowedError);
              CPPUNIT_ASSERT_DOUBLES_EQUAL(e3[i], extra.GetE3()[i], allowedError);
            }

            // This point should be just outside the 4 cube's outlet
            LatticeVector vLoc(2, 2, 5);
            // Create it
            VirtualSite vSite(initParams, extra, vLoc);

            // The iolet coords
            CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5, vSite.hv.posIolet.x, allowedError);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.5, vSite.hv.posIolet.y, allowedError);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.5, vSite.hv.posIolet.z, allowedError);

            // For D3Q15, there should be 5 cut links.
            CPPUNIT_ASSERT_EQUAL(size_t(5), vSite.neighbourGlobalIds.size());
            // each with q = 0.5
            CPPUNIT_ASSERT_DOUBLES_EQUAL(1.25, vSite.sumQiSq, allowedError);
            // For each site index in the neighbourhood
            for (Direction p = 0; p < lattice->GetNumVectors(); ++p)
            {
              //              if (lattice->GetVector(p).z < 0.)
              //              {
              //                CPPUNIT_ASSERT(vSite.streamingIndices[p] < 4 * 4 * 4 * 15);
              //              }
              //              else
              //              {
              //                CPPUNIT_ASSERT_EQUAL(site_t(4 * 4 * 4 * 15), vSite.streamingIndices[p]);
              //              }
            }

            // The velocity matrix = sum over i of:
            // x_i^2   x_i y_i x_i
            // x_i y_i y_i^2   y_i
            // x_i     y_i     1
            // distribn_t velMat[3][3] = { { 5.25, -1.25, 2.5 },
            //                             { -1.25, 5.25, -2.5 },
            //                             { 2.5, -2.5,  5.0 } };
            // It's inverse is
            distribn_t velMatInv[3][3] = { { 0.25, 0., -0.125 }, { 0., 0.25, 0.125 }, { -0.125,
                                                                                         0.125,
                                                                                         0.325 } };
            for (unsigned i = 0; i < 3; ++i)
            {
              for (unsigned j = 0; j < 3; ++j)
              {
                CPPUNIT_ASSERT_DOUBLES_EQUAL(velMatInv[i][j],
                                             vSite.velocityMatrixInv[i][j],
                                             allowedError);
              }
            }

            // Now check that appropriate entries have been added to the hydroCache
            CheckPointInCache(extra, LatticeVector(1, 1, 4), LatticePosition(1.5, -1.5, 0.5));
            CheckPointInCache(extra, LatticeVector(1, 3, 4), LatticePosition(1.5, 0.5, 0.5));
            CheckPointInCache(extra, LatticeVector(3, 1, 4), LatticePosition(-0.5, -1.5, 0.5));
            CheckPointInCache(extra, LatticeVector(3, 3, 4), LatticePosition(-0.5, 0.5, 0.5));
            CheckPointInCache(extra, LatticeVector(2, 2, 4), LatticePosition(0.5, -0.5, 0.5));
          }

          void TestStreamerInitialisation()
          {
            initParams.boundaryObject = outletBoundary;
            // Set up the ranges to cover Mid 3 (pure outlet) and Mid 5 (outlet/wall)
            initParams.siteRanges.resize(2);
            initParams.siteRanges[0].first = latDat->GetMidDomainCollisionCount(0)
                + latDat->GetMidDomainCollisionCount(1) + latDat->GetMidDomainCollisionCount(2);
            initParams.siteRanges[0].second = initParams.siteRanges[0].first
                + latDat->GetMidDomainCollisionCount(3);
            initParams.siteRanges[1].first = initParams.siteRanges[0].second
                + latDat->GetMidDomainCollisionCount(4);
            initParams.siteRanges[1].second = initParams.siteRanges[1].first
                + latDat->GetMidDomainCollisionCount(5);
            lb::streamers::VirtualSiteIolet<Collision> outletStreamer(initParams);

            // All the sites at the outlet plane (x, y, 3) should be in the cache.
            InOutLetCosine& outlet = *GetIolet(outletBoundary);
            VSExtra<Lattice>* extra = dynamic_cast<VSExtra<Lattice>*> (outlet.GetExtraData());
            CPPUNIT_ASSERT(extra != NULL);

            for (unsigned i = 1; i <= 4; ++i)
            {
              for (unsigned j = 1; j <= 4; ++j)
              {
                LatticeVector pos(i, j, 4);
                site_t globalIdx = latDat->GetGlobalNoncontiguousSiteIdFromGlobalCoords(pos);
                //                site_t localIdx = latDat->GetLocalContiguousIdFromGlobalNoncontiguousId(globalIdx);
                //                geometry::Site < geometry::LatticeData > site = latDat->GetSite(localIdx);

                RSHV::Map::iterator hvPtr = extra->hydroVarsCache.find(globalIdx);
                CPPUNIT_ASSERT(hvPtr != extra->hydroVarsCache.end());
              }
            }

            // And the reverse is true: every cache entry should be a site at the outlet plane
            for (RSHV::Map::iterator hvPtr = extra->hydroVarsCache.begin(); hvPtr
                != extra->hydroVarsCache.end(); ++hvPtr)
            {
              site_t globalIdx = hvPtr->first;
              LatticeVector pos;
              latDat->GetGlobalCoordsFromGlobalNoncontiguousSiteId(globalIdx, pos);
              CPPUNIT_ASSERT(hemelb::util::NumericalFunctions::IsInRange<LatticeCoordinate>(pos.x,
                                                                                            1,
                                                                                            4));
              CPPUNIT_ASSERT(hemelb::util::NumericalFunctions::IsInRange<LatticeCoordinate>(pos.y,
                                                                                            1,
                                                                                            4));
              CPPUNIT_ASSERT_EQUAL(LatticeCoordinate(4), pos.z);
            }
          }

          void TestStep()
          {
            initParams.boundaryObject = inletBoundary;
            lb::streamers::VirtualSiteIolet<Collision> inletStreamer(initParams);
            initParams.boundaryObject = outletBoundary;
            lb::streamers::VirtualSiteIolet<Collision> outletStreamer(initParams);

            InitialiseGradientHydroVars();

            // Stream and collide
            site_t offset = 0;
            offset += latDat->GetMidDomainCollisionCount(0);
            offset += latDat->GetMidDomainCollisionCount(1);
            inletStreamer.DoStreamAndCollide<false> (offset,
                                                     latDat->GetMidDomainCollisionCount(2),
                                                     lbmParams,
                                                     static_cast<geometry::LatticeData*> (latDat),
                                                     *propertyCache);
            offset += latDat->GetMidDomainCollisionCount(2);

            outletStreamer.StreamAndCollide<false> (offset,
                                                    latDat->GetMidDomainCollisionCount(3),
                                                    lbmParams,
                                                    latDat,
                                                    *propertyCache);
            offset += latDat->GetMidDomainCollisionCount(3);

            inletStreamer.StreamAndCollide<false> (offset,
                                                   latDat->GetMidDomainCollisionCount(4),
                                                   lbmParams,
                                                   latDat,
                                                   *propertyCache);
            offset += latDat->GetMidDomainCollisionCount(4);

            outletStreamer.StreamAndCollide<false> (offset,
                                                    latDat->GetMidDomainCollisionCount(5),
                                                    lbmParams,
                                                    latDat,
                                                    *propertyCache);

            // Now every entry in the RSHV cache should have been updated
            CheckAllHVUpdated(inletBoundary, 1);
            CheckAllHVUpdated(outletBoundary, 1);

            // Stream and collide
            offset = 0;
            offset += latDat->GetMidDomainCollisionCount(0);
            offset += latDat->GetMidDomainCollisionCount(1);
            inletStreamer.DoPostStep<false> (offset,
                                             latDat->GetMidDomainCollisionCount(2),
                                             lbmParams,
                                             static_cast<geometry::LatticeData*> (latDat),
                                             *propertyCache);
            offset += latDat->GetMidDomainCollisionCount(2);

            outletStreamer.DoPostStep<false> (offset,
                                              latDat->GetMidDomainCollisionCount(3),
                                              lbmParams,
                                              latDat,
                                              *propertyCache);
            offset += latDat->GetMidDomainCollisionCount(3);

            inletStreamer.DoPostStep<false> (offset,
                                             latDat->GetMidDomainCollisionCount(4),
                                             lbmParams,
                                             latDat,
                                             *propertyCache);
            offset += latDat->GetMidDomainCollisionCount(4);

            outletStreamer.DoPostStep<false> (offset,
                                              latDat->GetMidDomainCollisionCount(5),
                                              lbmParams,
                                              latDat,
                                              *propertyCache);

            // Check that all the vsites have sensible hydro values
            InOutLetCosine* inlet = GetIolet(inletBoundary);
            VSExtra<Lattice> * inExtra = dynamic_cast<VSExtra<Lattice>*> (inlet->GetExtraData());

            for (VirtualSite::Map::iterator vsIt = inExtra->vSites.begin(); vsIt
                != inExtra->vSites.end(); ++vsIt)
            {
              site_t vSiteGlobalIdx = vsIt->first;
              VirtualSite& vSite = vsIt->second;

              CPPUNIT_ASSERT_EQUAL(LatticeTimeStep(1), vSite.hv.t);
              CPPUNIT_ASSERT_DOUBLES_EQUAL(LatticeDensity(1.045), vSite.hv.rho, allowedError);

              LatticeVector pos;
              latDat->GetGlobalCoordsFromGlobalNoncontiguousSiteId(vSiteGlobalIdx, pos);

              if (vSite.neighbourGlobalIds.size() > 3)
                CPPUNIT_ASSERT_DOUBLES_EQUAL(GetVelocity(pos).z, vSite.hv.u.z, allowedError);

            }

            InOutLetCosine* outlet = GetIolet(outletBoundary);
            VSExtra<Lattice> * outExtra = dynamic_cast<VSExtra<Lattice>*> (outlet->GetExtraData());

            for (VirtualSite::Map::iterator vsIt = outExtra->vSites.begin(); vsIt
                != outExtra->vSites.end(); ++vsIt)
            {
              site_t vSiteGlobalIdx = vsIt->first;
              VirtualSite& vSite = vsIt->second;

              CPPUNIT_ASSERT_EQUAL(LatticeTimeStep(1), vSite.hv.t);
              CPPUNIT_ASSERT_DOUBLES_EQUAL(LatticeDensity(0.995), vSite.hv.rho, allowedError);

              LatticeVector pos;
              latDat->GetGlobalCoordsFromGlobalNoncontiguousSiteId(vSiteGlobalIdx, pos);

              if (vSite.neighbourGlobalIds.size() > 3)
                CPPUNIT_ASSERT_DOUBLES_EQUAL(GetVelocity(pos).z, vSite.hv.u.z, allowedError);

            }
          }
        private:
          lb::MacroscopicPropertyCache* propertyCache;
          lb::iolets::BoundaryValues* outletBoundary;
          lb::iolets::BoundaryValues* inletBoundary;
          lb::lattices::LatticeInfo* lattice;

          void CheckAllHVUpdated(lb::iolets::BoundaryValues* iolets, LatticeTimeStep expectedT)
          {
            VSExtra<Lattice> * extra =
                dynamic_cast<VSExtra<Lattice>*> (iolets->GetLocalIolet(0)->GetExtraData());
            for (RSHV::Map::iterator hvPtr = extra->hydroVarsCache.begin(); hvPtr
                != extra->hydroVarsCache.end(); ++hvPtr)
            {
              site_t siteGlobalIdx = hvPtr->first;
              LatticeVector sitePos;
              latDat->GetGlobalCoordsFromGlobalNoncontiguousSiteId(siteGlobalIdx, sitePos);
              RSHV& hv = hvPtr->second;
              CPPUNIT_ASSERT_EQUAL(expectedT, hv.t);
              CPPUNIT_ASSERT_DOUBLES_EQUAL(GetDensity(sitePos), hv.rho, allowedError);
              LatticeVelocity u = GetVelocity(sitePos);
              for (unsigned i = 0; i < 3; ++i)
                CPPUNIT_ASSERT_DOUBLES_EQUAL(u[i], hv.u[i], allowedError);

            }
          }

          void CheckPointInCache(LocalVSExtra& extra, LatticeVector expectedPt,
                                 LatticePosition expectedIoletPos)
          {
            site_t expectedGlobalIdx =
                latDat->GetGlobalNoncontiguousSiteIdFromGlobalCoords(expectedPt);
            RSHV::Map::iterator hvPtr = extra.hydroVarsCache.find(expectedGlobalIdx);

            CPPUNIT_ASSERT(hvPtr != extra.hydroVarsCache.end());
            for (unsigned i = 0; i < 3; ++i)
              CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedIoletPos[i],
                                           hvPtr->second.posIolet[i],
                                           allowedError);

          }

          static LatticeVelocity GetVelocity(LatticeVector pos)
          {
            LatticeVelocity c(0.01, 0.01, 0.);
            LatticeVelocity u(0.);
            u.z = c.Dot(pos);
            return u;
          }
          static LatticeDensity GetDensity(LatticeVector pos)
          {
            util::Vector3D<LatticeDensity> gradRho(0., 0., -1e-2);
            return 1.045 + gradRho.Dot(pos);
          }
          void InitialiseGradientHydroVars()
          {
            for (unsigned i = 1; i <= 4; ++i)
            {
              for (unsigned j = 1; j <= 4; ++j)
              {
                for (unsigned k = 1; k <= 4; ++k)
                {
                  LatticeVector pos(i, j, k);
                  site_t siteIdx = latDat->GetContiguousSiteId(pos);
                  //geometry::Site < geometry::LatticeData > site = latDat->GetSite(siteIdx);
                  distribn_t* fOld = latDat->GetFNew(siteIdx * Lattice::NUMVECTORS);
                  LatticeDensity rho = GetDensity(pos);
                  LatticeVelocity u = GetVelocity(pos);
                  u *= rho;
                  Lattice::CalculateFeq(rho, u.x, u.y, u.z, fOld);
                }
              }
            }
            latDat->SwapOldAndNew();
          }
      };
      CPPUNIT_TEST_SUITE_REGISTRATION( VirtualSiteIoletStreamerTests);
    }
  }
}

#endif
