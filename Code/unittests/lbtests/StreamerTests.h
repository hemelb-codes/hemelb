
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_LBTESTS_STREAMERTESTS_H
#define HEMELB_UNITTESTS_LBTESTS_STREAMERTESTS_H

#include <cppunit/TestFixture.h>
#include <iostream>
#include <sstream>

#include "lb/streamers/Streamers.h"
#include "geometry/SiteData.h"

#include "unittests/helpers/FourCubeBasedTestFixture.h"

namespace hemelb
{
  namespace unittests
  {
    namespace lbtests
    {
      static const distribn_t allowedError = 1e-10;

      /**
       * StreamerTests:
       *
       * This class tests the streamer implementations. We assume the collision operators are
       * correct (as they're tested elsewhere), then compare the post-streamed values with
       * the values we expect to have been streamed there.
       */
      class StreamerTests : public helpers::FourCubeBasedTestFixture
      {
          CPPUNIT_TEST_SUITE ( StreamerTests);
          CPPUNIT_TEST ( TestSimpleCollideAndStream);
          CPPUNIT_TEST ( TestBouzidiFirdaousLallemand);
          CPPUNIT_TEST ( TestSimpleBounceBack);
          CPPUNIT_TEST ( TestGuoZhengShi);
          CPPUNIT_TEST ( TestNashZerothOrderPressureIolet);
          CPPUNIT_TEST ( TestNashZerothOrderPressureBB);
          CPPUNIT_TEST ( TestJunkYangEquivalentToBounceBack);CPPUNIT_TEST_SUITE_END();
        public:

          void setUp()
          {

            FourCubeBasedTestFixture::setUp();
            propertyCache = new lb::MacroscopicPropertyCache(*simState, *latDat);
            normalCollision
                = new lb::collisions::Normal<lb::kernels::LBGK<lb::lattices::D3Q15> >(initParams);
          }

          void tearDown()
          {
            delete propertyCache;
            delete normalCollision;

            FourCubeBasedTestFixture::tearDown();
          }

          void TestSimpleCollideAndStream()
          {
            lb::streamers::SimpleCollideAndStream<lb::collisions::Normal<lb::kernels::LBGK<
                lb::lattices::D3Q15> > > simpleCollideAndStream(initParams);

            // Initialise fOld in the lattice data. We choose values so that each site has
            // an anisotropic distribution function, and that each site's function is
            // distinguishable.
            LbTestsHelper::InitialiseAnisotropicTestData<lb::lattices::D3Q15>(latDat);

            // Use the streaming operator on the entire lattice.
            simpleCollideAndStream.StreamAndCollide<false> (0,
                                                            latDat->GetLocalFluidSiteCount(),
                                                            lbmParams,
                                                            latDat,
                                                            *propertyCache);

            // Now, go over each lattice site and check each value in f_new is correct.
            for (site_t streamedToSite = 0; streamedToSite < latDat->GetLocalFluidSiteCount(); ++streamedToSite)
            {
              geometry::Site < geometry::LatticeData > streamedSite
                  = latDat->GetSite(streamedToSite);

              distribn_t* streamedToFNew = latDat->GetFNew(lb::lattices::D3Q15::NUMVECTORS
                  * streamedToSite);

              for (unsigned int streamedDirection = 0; streamedDirection
                  < lb::lattices::D3Q15::NUMVECTORS; ++streamedDirection)
              {

                site_t
                    streamerIndex =
                        streamedSite.GetStreamedIndex<lb::lattices::D3Q15> (lb::lattices::D3Q15::INVERSEDIRECTIONS[streamedDirection]);

                // If this site streamed somewhere sensible, it must have been streamed to.
                if (streamerIndex >= 0 && streamerIndex < (lb::lattices::D3Q15::NUMVECTORS
                    * latDat->GetLocalFluidSiteCount()))
                {
                  site_t streamerSiteId = streamerIndex / lb::lattices::D3Q15::NUMVECTORS;

                  // Calculate streamerFOld at this site.
                  distribn_t streamerFOld[lb::lattices::D3Q15::NUMVECTORS];
                  LbTestsHelper::InitialiseAnisotropicTestData<lb::lattices::D3Q15>(streamerSiteId,
                                                                                    streamerFOld);

                  // Calculate what the value streamed to site streamedToSite should be.
                  lb::kernels::HydroVars<lb::kernels::LBGK<lb::lattices::D3Q15> >
                      streamerHydroVars(streamerFOld);
                  normalCollision->CalculatePreCollision(streamerHydroVars, streamedSite);

                  normalCollision->Collide(lbmParams, streamerHydroVars);

                  // F_new should be equal to the value that was streamed from this other site
                  // in the same direction as we're streaming from.
                  CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("SimpleCollideAndStream, StreamAndCollide",
                                                       streamerHydroVars.GetFPostCollision()[streamedDirection],
                                                       streamedToFNew[streamedDirection],
                                                       allowedError);
                }
              }
            }
          }

          void TestBouzidiFirdaousLallemand()
          {
            // Initialise fOld in the lattice data. We choose values so that each site has
            // an anisotropic distribution function, and that each site's function is
            // distinguishable.
            LbTestsHelper::InitialiseAnisotropicTestData<lb::lattices::D3Q15>(latDat);
            lb::streamers::BouzidiFirdaousLallemand<lb::collisions::Normal<lb::kernels::LBGK<
                lb::lattices::D3Q15> > >::Type bfl(initParams);

            bfl.StreamAndCollide<false> (0,
                                         latDat->GetLocalFluidSiteCount(),
                                         lbmParams,
                                         latDat,
                                         *propertyCache);
            bfl.PostStep<false> (0,
                                 latDat->GetLocalFluidSiteCount(),
                                 lbmParams,
                                 latDat,
                                 *propertyCache);

            // Now, go over each lattice site and check each value in f_new is correct.
            for (site_t streamedToSite = 0; streamedToSite < latDat->GetLocalFluidSiteCount(); ++streamedToSite)
            {
              const geometry::Site<geometry::LatticeData> streamedSite =
                  latDat->GetSite(streamedToSite);

              distribn_t* streamedToFNew = latDat->GetFNew(lb::lattices::D3Q15::NUMVECTORS
                  * streamedToSite);

              for (unsigned int streamedDirection = 0; streamedDirection
                  < lb::lattices::D3Q15::NUMVECTORS; ++streamedDirection)
              {
                unsigned int oppDirection =
                    lb::lattices::D3Q15::INVERSEDIRECTIONS[streamedDirection];

                site_t streamerIndex =
                    streamedSite.GetStreamedIndex<lb::lattices::D3Q15> (oppDirection);

                geometry::Site < geometry::LatticeData > streamerSite
                    = latDat->GetSite(streamerIndex);

                // If this site streamed somewhere sensible, it must have been streamed to.
                if (streamerIndex >= 0 && streamerIndex < (lb::lattices::D3Q15::NUMVECTORS
                    * latDat->GetLocalFluidSiteCount()))
                {
                  site_t streamerSiteId = streamerIndex / lb::lattices::D3Q15::NUMVECTORS;

                  // Calculate streamerFOld at this site.
                  distribn_t streamerFOld[lb::lattices::D3Q15::NUMVECTORS];
                  LbTestsHelper::InitialiseAnisotropicTestData<lb::lattices::D3Q15>(streamerSiteId,
                                                                                    streamerFOld);

                  // Calculate what the value streamed to site streamedToSite should be.
                  lb::kernels::HydroVars<lb::kernels::LBGK<lb::lattices::D3Q15> >
                      streamerHydroVars(streamerFOld);
                  normalCollision->CalculatePreCollision(streamerHydroVars, streamerSite);

                  normalCollision->Collide(lbmParams, streamerHydroVars);

                  // F_new should be equal to the value that was streamed from this other site
                  // in the same direction as we're streaming from.
                  CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("SimpleCollideAndStream, StreamAndCollide",
                                                       streamerHydroVars.GetFPostCollision()[streamedDirection],
                                                       streamedToFNew[streamedDirection],
                                                       allowedError);
                }
                else if (streamedSite.GetSiteType() == geometry::INLET_TYPE
                    || streamedSite.GetSiteType() == geometry::OUTLET_TYPE)
                {
                  // No reason to further test an inlet/outlet site.
                  // Pass.
                }
                else
                {
                  std::stringstream message;
                  message << "Site: " << streamedToSite << " Direction " << oppDirection
                      << " Data: " << streamedSite.GetSiteData().GetWallIntersectionData()
                      << std::flush;
                  CPPUNIT_ASSERT_MESSAGE("Expected to find a boundary "
                                           "opposite an unstreamed-to direction " + message.str(),
                                         streamedSite.HasWall(oppDirection));
                  // Test disabled due to RegressionTests issue, see discussion in #87
                  CPPUNIT_ASSERT_MESSAGE("Expect defined cut distance opposite an unstreamed-to direction "
                                             + message.str(),
                                         streamedSite.GetWallDistance<lb::lattices::D3Q15> (oppDirection)
                                             != -1.0);

                  // To verify the operation of the BFL boundary condition, we'll need:
                  // - the distance to the wall * 2
                  distribn_t twoQ = 2.0
                      * streamedSite.GetWallDistance<lb::lattices::D3Q15> (oppDirection);

                  // - the post-collision distribution at the current site.
                  distribn_t streamedToSiteFOld[lb::lattices::D3Q15::NUMVECTORS];

                  // (initialise it to f_old).
                  LbTestsHelper::InitialiseAnisotropicTestData<lb::lattices::D3Q15>(streamedToSite,
                                                                                    streamedToSiteFOld);

                  lb::kernels::HydroVars<lb::kernels::LBGK<lb::lattices::D3Q15> >
                      hydroVars(streamedToSiteFOld);

                  normalCollision->CalculatePreCollision(hydroVars, streamedSite);

                  // (find post-collision values using the collision operator).
                  normalCollision->Collide(lbmParams, hydroVars);

                  // - finally, the post-collision distribution at the site which is one further
                  // away from the wall in this direction.
                  distribn_t awayFromWallFOld[lb::lattices::D3Q15::NUMVECTORS];

                  site_t awayFromWallIndex =
                      streamedSite.GetStreamedIndex<lb::lattices::D3Q15> (streamedDirection)
                          / lb::lattices::D3Q15::NUMVECTORS;

                  // If there's a valid index in that direction, use BFL
                  if (awayFromWallIndex >= 0 && awayFromWallIndex
                      < latDat->GetLocalFluidSiteCount())
                  {
                    const geometry::Site<geometry::LatticeData> awayFromWallSite =
                        latDat->GetSite(awayFromWallIndex);

                    // (initialise it to f_old).
                    LbTestsHelper::InitialiseAnisotropicTestData<lb::lattices::D3Q15>(awayFromWallIndex,
                                                                                      awayFromWallFOld);

                    lb::kernels::HydroVars<lb::kernels::LBGK<lb::lattices::D3Q15> >
                        awayFromWallsHydroVars(awayFromWallFOld);

                    normalCollision->CalculatePreCollision(awayFromWallsHydroVars, awayFromWallSite);

                    // (find post-collision values using the collision operator).
                    normalCollision->Collide(lbmParams, awayFromWallsHydroVars);

                    distribn_t toWallOld = hydroVars.GetFPostCollision()[oppDirection];

                    distribn_t toWallNew = awayFromWallsHydroVars.GetFPostCollision()[oppDirection];

                    distribn_t oppWallOld = hydroVars.GetFPostCollision()[streamedDirection];

                    // The streamed value should be as given below.
                    distribn_t streamed = (twoQ < 1.0)
                      ? (toWallNew + twoQ * (toWallOld - toWallNew))
                      : (oppWallOld + (1. / twoQ) * (toWallOld - oppWallOld));

                    std::stringstream msg(std::stringstream::in);
                    msg << "BouzidiFirdaousLallemand, PostStep: site " << streamedToSite
                        << " direction " << streamedDirection;

                    // Assert that this is the case.
                    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(msg.str(),
                                                         streamed,
                                                         * (latDat->GetFNew(streamedToSite
                                                             * lb::lattices::D3Q15::NUMVECTORS
                                                             + streamedDirection)),
                                                         allowedError);
                  }

                  // With no valid lattice site, simple bounce-back will be performed.
                  else
                  {
                    std::stringstream msg(std::stringstream::in);
                    msg << "BouzidiFirdaousLallemand, PostStep by simple bounce-back: site "
                        << streamedToSite << " direction " << streamedDirection;

                    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(msg.str(),
                                                         hydroVars.GetFPostCollision()[oppDirection],
                                                         * (latDat->GetFNew(streamedToSite
                                                             * lb::lattices::D3Q15::NUMVECTORS
                                                             + streamedDirection)),
                                                         allowedError);
                  }
                }
              }
            }
          }

          void TestSimpleBounceBack()
          {
            // Initialise fOld in the lattice data. We choose values so that each site has
            // an anisotropic distribution function, and that each site's function is
            // distinguishable.
            LbTestsHelper::InitialiseAnisotropicTestData<lb::lattices::D3Q15>(latDat);

            site_t firstWallSite = latDat->GetMidDomainCollisionCount(0);
            site_t wallSitesCount = latDat->GetMidDomainCollisionCount(1);

            // Check that the lattice has the expected number of sites labeled as pure wall (otherwise this test is void)
            CPPUNIT_ASSERT_EQUAL_MESSAGE("Unexpected number of pure wall sites.",
                                         site_t(24),
                                         latDat->GetMidDomainCollisionCount(1));

            site_t offset = 0;

            // Mid-Fluid sites use simple collide and stream
            lb::streamers::SimpleCollideAndStream<lb::collisions::Normal<lb::kernels::LBGK<
                lb::lattices::D3Q15> > > simpleCollideAndStream(initParams);

            simpleCollideAndStream.StreamAndCollide<false> (offset,
                                                            latDat->GetMidDomainCollisionCount(0),
                                                            lbmParams,
                                                            latDat,
                                                            *propertyCache);
            offset += latDat->GetMidDomainCollisionCount(0);

            // Wall sites use simple bounce back
            lb::streamers::SimpleBounceBack<lb::collisions::Normal<lb::kernels::LBGK<
                lb::lattices::D3Q15> > >::Type simpleBounceBack(initParams);

            simpleBounceBack.StreamAndCollide<false> (offset,
                                                      latDat->GetMidDomainCollisionCount(1),
                                                      lbmParams,
                                                      latDat,
                                                      *propertyCache);
            offset += latDat->GetMidDomainCollisionCount(1);

            // Consider inlet/outlets and their walls as mid-fluid sites
            simpleCollideAndStream.StreamAndCollide<false> (offset,
                                                            latDat->GetLocalFluidSiteCount()
                                                                - offset,
                                                            lbmParams,
                                                            latDat,
                                                            *propertyCache);
            offset += latDat->GetLocalFluidSiteCount() - offset;

            // Sanity check
            CPPUNIT_ASSERT_EQUAL_MESSAGE("Total number of sites",
                                         offset,
                                         latDat->GetLocalFluidSiteCount());

            /*
             *  Loop over the wall sites and check whether they got properly streamed on or bounced back
             *  depending on where they sit relative to the wall. We ignore mid-Fluid sites since
             *  StreamAndCollide was tested before.
             */
            for (site_t wallSiteLocalIndex = 0; wallSiteLocalIndex < wallSitesCount; wallSiteLocalIndex++)
            {
              site_t streamedToSite = firstWallSite + wallSiteLocalIndex;
              const geometry::Site<geometry::LatticeData> streamedSite =
                  latDat->GetSite(streamedToSite);
              distribn_t* streamedToFNew = latDat->GetFNew(lb::lattices::D3Q15::NUMVECTORS
                  * streamedToSite);

              for (unsigned int streamedDirection = 0; streamedDirection
                  < lb::lattices::D3Q15::NUMVECTORS; ++streamedDirection)
              {
                unsigned oppDirection = lb::lattices::D3Q15::INVERSEDIRECTIONS[streamedDirection];

                // Index of the site streaming to streamedToSite via direction streamedDirection
                site_t streamerIndex =
                    streamedSite.GetStreamedIndex<lb::lattices::D3Q15> (oppDirection);

                // Is streamerIndex a valid index?
                if (streamerIndex >= 0 && streamerIndex < (lb::lattices::D3Q15::NUMVECTORS
                    * latDat->GetLocalFluidSiteCount()))
                {
                  // The streamer index is a valid index in the domain, therefore stream and collide has happened
                  site_t streamerSiteId = streamerIndex / lb::lattices::D3Q15::NUMVECTORS;

                  // Calculate streamerFOld at this site.
                  distribn_t streamerFOld[lb::lattices::D3Q15::NUMVECTORS];
                  LbTestsHelper::InitialiseAnisotropicTestData<lb::lattices::D3Q15>(streamerSiteId,
                                                                                    streamerFOld);

                  // Calculate what the value streamed to site streamedToSite should be.
                  lb::kernels::HydroVars<lb::kernels::LBGK<lb::lattices::D3Q15> >
                      streamerHydroVars(streamerFOld);
                  normalCollision->CalculatePreCollision(streamerHydroVars, streamedSite);

                  normalCollision->Collide(lbmParams, streamerHydroVars);

                  // F_new should be equal to the value that was streamed from this other site
                  // in the same direction as we're streaming from.
                  CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("SimpleCollideAndStream, StreamAndCollide",
                                                       streamerHydroVars.GetFPostCollision()[streamedDirection],
                                                       streamedToFNew[streamedDirection],
                                                       allowedError);
                }
                else
                {
                  // The streamer index shows that no one has streamed to streamedToSite direction
                  // streamedDirection, therefore bounce back has happened in that site for that direction

                  // Initialise streamedToSiteFOld with the original data
                  distribn_t streamerToSiteFOld[lb::lattices::D3Q15::NUMVECTORS];
                  LbTestsHelper::InitialiseAnisotropicTestData<lb::lattices::D3Q15>(streamedToSite,
                                                                                    streamerToSiteFOld);
                  lb::kernels::HydroVars<lb::kernels::LBGK<lb::lattices::D3Q15> >
                      hydroVars(streamerToSiteFOld);
                  normalCollision->CalculatePreCollision(hydroVars, streamedSite);

                  // Simulate post-collision using the collision operator.
                  normalCollision->Collide(lbmParams, hydroVars);

                  // After streaming FNew in a given direction must be f post-collision in the opposite direction
                  // following collision
                  std::stringstream msg(std::stringstream::in);
                  msg << "Simple bounce-back: site " << streamedToSite << " direction "
                      << streamedDirection;
                  CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(msg.str(),
                                                       streamedToFNew[streamedDirection],
                                                       hydroVars.GetFPostCollision()[oppDirection],
                                                       allowedError);
                }
              }
            }
          }

          void TestGuoZhengShi()
          {
            lb::streamers::GuoZhengShi<lb::collisions::Normal<
                lb::kernels::LBGK<lb::lattices::D3Q15> > >::Type guoZhengShi(initParams);

            for (double assignedWallDistance = 0.4; assignedWallDistance < 1.0; assignedWallDistance
                += 0.5)
            {
              // Initialise fOld in the lattice data. We choose values so that each site has
              // an anisotropic distribution function, and that each site's function is
              // distinguishable.
              LbTestsHelper::InitialiseAnisotropicTestData<lb::lattices::D3Q15>(latDat);

              // Make some fairly arbitrary choices early on.
              const site_t chosenSite = 0;
              const geometry::Site<geometry::LatticeData>& streamer = latDat->GetSite(chosenSite);

              const Direction chosenUnstreamedDirection = 5;
              const Direction chosenWallDirection =
                  lb::lattices::D3Q15::INVERSEDIRECTIONS[chosenUnstreamedDirection];
              const Direction chosenDoubleWallDirection1 = 7;
              const Direction chosenDoubleWallDirection2 = 8;

              // Enforce that there's a boundary in the wall direction.
              latDat->SetHasWall(chosenSite, chosenWallDirection);
              latDat->SetHasWall(chosenSite, chosenDoubleWallDirection1);
              latDat->SetHasWall(chosenSite, chosenDoubleWallDirection2);
              latDat->SetBoundaryDistance(chosenSite, chosenWallDirection, assignedWallDistance);
              latDat->SetBoundaryDistance(chosenSite,
                                          chosenDoubleWallDirection1,
                                          assignedWallDistance);
              latDat->SetBoundaryDistance(chosenSite,
                                          chosenDoubleWallDirection2,
                                          assignedWallDistance);

              // Perform the collision and streaming.
              guoZhengShi.StreamAndCollide<false> (chosenSite, 1, lbmParams, latDat, *propertyCache);

              // Calculate the distributions at the chosen site up to post-collision.
              distribn_t streamerFOld[lb::lattices::D3Q15::NUMVECTORS];
              LbTestsHelper::InitialiseAnisotropicTestData<lb::lattices::D3Q15>(chosenSite,
                                                                                streamerFOld);
              lb::kernels::HydroVars<lb::kernels::LBGK<lb::lattices::D3Q15> >
                  streamerHydroVars(streamerFOld);
              normalCollision->CalculatePreCollision(streamerHydroVars, streamer);
              normalCollision->Collide(lbmParams, streamerHydroVars);

              // Check each streamed direction.
              for (Direction streamedDirection = 0; streamedDirection
                  < lb::lattices::D3Q15::NUMVECTORS; ++streamedDirection)
              {
                switch (streamedDirection)
                {
                  case 5:
                  {
                    // We point away from the wall and there is a fluid site out in this direction
                    // So we have regular GZS.

                    LatticeVelocity velocityWall;
                    distribn_t fNeqWall;

                    // This is the first means of estimating from the source paper: only
                    // use the nearest fluid site.
                    LatticeVelocity velocityEstimate1 = streamerHydroVars.momentum * (1. - 1.
                        / assignedWallDistance) / streamerHydroVars.density;

                    distribn_t fNeqEstimate1 = streamerHydroVars.GetFNeq()[streamedDirection];

                    if (assignedWallDistance < 0.75)
                    {
                      // This is the second method for estimating: using the next fluid site
                      // away from the wall.
                      const site_t nextSiteAwayFromWall = streamer.GetStreamedIndex<
                          lb::lattices::D3Q15> (streamedDirection)
                          / lb::lattices::D3Q15::NUMVECTORS;
                      const geometry::Site<geometry::LatticeData>& nextSiteAway =
                          latDat->GetSite(nextSiteAwayFromWall);
                      distribn_t nextSiteOutFOld[lb::lattices::D3Q15::NUMVECTORS];
                      LbTestsHelper::InitialiseAnisotropicTestData<lb::lattices::D3Q15>(nextSiteAwayFromWall,
                                                                                        nextSiteOutFOld);
                      lb::kernels::HydroVars<lb::kernels::LBGK<lb::lattices::D3Q15> >
                          nextSiteOutHydroVars(nextSiteOutFOld);
                      normalCollision->CalculatePreCollision(nextSiteOutHydroVars, nextSiteAway);

                      LatticeVelocity velocityEstimate2 = nextSiteOutHydroVars.momentum
                          * (assignedWallDistance - 1.) / ( (1. + assignedWallDistance)
                          * nextSiteOutHydroVars.density);

                      distribn_t fNeqEstimate2 = nextSiteOutHydroVars.GetFNeq()[streamedDirection];

                      // The actual value is taken to be an interpolation between the two
                      // estimates.
                      velocityWall = velocityEstimate1 * assignedWallDistance + velocityEstimate2
                          * (1. - assignedWallDistance);

                      fNeqWall = assignedWallDistance * fNeqEstimate1 + (1. - assignedWallDistance)
                          * fNeqEstimate2;
                    }
                    else
                    {
                      // Only use the first estimate
                      velocityWall = velocityEstimate1;
                      fNeqWall = fNeqEstimate1;
                    }

                    util::Vector3D<distribn_t> momentumWall = velocityWall
                        * streamerHydroVars.density;

                    // Compute the wall equilibrium dist
                    distribn_t fEqm[lb::lattices::D3Q15::NUMVECTORS];
                    lb::lattices::D3Q15::CalculateFeq(streamerHydroVars.density,
                                                      momentumWall.x,
                                                      momentumWall.y,
                                                      momentumWall.z,
                                                      fEqm);
                    // Perform collision on the wall f's
                    distribn_t prediction = fEqm[streamedDirection] + (1.0 + lbmParams->GetOmega())
                        * fNeqWall;
                    // This is the answer from the code we're testing
                    distribn_t streamedFNew = latDat->GetFNew(lb::lattices::D3Q15::NUMVECTORS
                        * chosenSite)[streamedDirection];

                    CPPUNIT_ASSERT_DOUBLES_EQUAL(prediction, streamedFNew, allowedError);
                    break;
                  }
                  case 7:
                  case 8:
                    // These guys both point at the wall and are opposite each other.
                    // We might do SBB or GZS depending on how close the wall is.
                    if (assignedWallDistance < 0.75)
                    {
                      // It's SBB
                      Direction inv = lb::lattices::D3Q15::INVERSEDIRECTIONS[streamedDirection];
                      distribn_t prediction = streamerHydroVars.GetFPostCollision()[inv];
                      // This is the answer from the code we're testing
                      distribn_t streamedFNew = latDat->GetFNew(lb::lattices::D3Q15::NUMVECTORS
                          * chosenSite)[streamedDirection];
                      CPPUNIT_ASSERT_DOUBLES_EQUAL(prediction, streamedFNew, allowedError);
                    }
                    else
                    {
                      // It's GZS with extrapolation from this site only

                      // This is the first means of estimating from the source paper: only
                      // use the nearest fluid site.
                      LatticeVelocity velocityWall = streamerHydroVars.momentum * (1. - 1.
                          / assignedWallDistance) / streamerHydroVars.density;

                      distribn_t fNeqWall = streamerHydroVars.GetFNeq()[streamedDirection];

                      util::Vector3D<distribn_t> momentumWall = velocityWall
                          * streamerHydroVars.density;

                      // Compute the wall equilibrium dist
                      distribn_t fEqm[lb::lattices::D3Q15::NUMVECTORS];
                      lb::lattices::D3Q15::CalculateFeq(streamerHydroVars.density,
                                                        momentumWall.x,
                                                        momentumWall.y,
                                                        momentumWall.z,
                                                        fEqm);

                      // Perform collision on the wall f's
                      distribn_t prediction = fEqm[streamedDirection] + (1.0
                          + lbmParams->GetOmega()) * fNeqWall;
                      // This is the answer from the code we're testing
                      distribn_t streamedFNew = latDat->GetFNew(lb::lattices::D3Q15::NUMVECTORS
                          * chosenSite)[streamedDirection];

                      CPPUNIT_ASSERT_DOUBLES_EQUAL(prediction, streamedFNew, allowedError);

                    }
                    break;
                  case 6:
                    // We point towards the wall.
                    break;

                  default:
                    // We have nothing to do with a wall so simple streaming
                    const site_t streamedIndex =
                        streamer.GetStreamedIndex<lb::lattices::D3Q15> (streamedDirection);
                    distribn_t streamedToFNew = *latDat->GetFNew(streamedIndex);

                    // F_new should be equal to the value that was streamed from this other site
                    // in the same direction as we're streaming from.
                    CPPUNIT_ASSERT_DOUBLES_EQUAL(streamerHydroVars.GetFPostCollision()[streamedDirection],
                                                 streamedToFNew,
                                                 allowedError);
                    break;
                }

              }
            }
          }

          /**
           * Junk&Yang should behave like simple bounce back when fluid sites are 0.5 lattice length units away
           * from the domain boundary.
           */
          void TestJunkYangEquivalentToBounceBack()
          {
            // Initialise fOld in the lattice data. We choose values so that each site has
            // an anisotropic distribution function, and that each site's function is
            // distinguishable.
            LbTestsHelper::InitialiseAnisotropicTestData<lb::lattices::D3Q15>(latDat);
            // Setting all the wall distances to 0.5 will make Junk&Yang behave like Simple Bounce Back
            LbTestsHelper::SetWallAndIoletDistances<lb::lattices::D3Q15>(*latDat, 0.5);

            site_t firstWallSite = latDat->GetMidDomainCollisionCount(0);
            site_t wallSitesCount = latDat->GetMidDomainCollisionCount(1) - firstWallSite;

            // Check that the lattice has sites labeled as wall (otherwise this test is void)
            CPPUNIT_ASSERT_EQUAL_MESSAGE("Number of non-wall sites", wallSitesCount, (site_t) 16);

            site_t offset = 0;

            // Mid-Fluid sites use simple collide and stream
            lb::streamers::SimpleCollideAndStream<lb::collisions::Normal<lb::kernels::LBGK<
                lb::lattices::D3Q15> > > simpleCollideAndStream(initParams);

            simpleCollideAndStream.StreamAndCollide<false> (offset,
                                                            latDat->GetMidDomainCollisionCount(0),
                                                            lbmParams,
                                                            latDat,
                                                            *propertyCache);
            offset += latDat->GetMidDomainCollisionCount(0);

            // Wall sites use junk and yang
            initParams.siteRanges.push_back(std::pair<site_t, site_t>(offset,
                                                                      offset
                                                                          + latDat->GetMidDomainCollisionCount(1)));
            lb::streamers::JunkYang<lb::collisions::Normal<lb::kernels::LBGK<lb::lattices::D3Q15> > >::Type
                junkYang(initParams);

            junkYang.StreamAndCollide<false> (offset,
                                              latDat->GetMidDomainCollisionCount(1),
                                              lbmParams,
                                              latDat,
                                              *propertyCache);

            junkYang.PostStep<false> (offset,
                                      latDat->GetMidDomainCollisionCount(1),
                                      lbmParams,
                                      latDat,
                                      *propertyCache);

            offset += latDat->GetMidDomainCollisionCount(1);

            // Consider inlet/outlets and their walls as mid-fluid sites
            simpleCollideAndStream.StreamAndCollide<false> (offset,
                                                            latDat->GetLocalFluidSiteCount()
                                                                - offset,
                                                            lbmParams,
                                                            latDat,
                                                            *propertyCache);
            offset += latDat->GetLocalFluidSiteCount() - offset;

            // Sanity check
            CPPUNIT_ASSERT_EQUAL_MESSAGE("Total number of sites",
                                         offset,
                                         latDat->GetLocalFluidSiteCount());

            /*
             *  Loop over the wall sites and check whether they got properly streamed on or bounced back
             *  depending on where they sit relative to the wall. We ignore mid-Fluid sites since
             *  StreamAndCollide was tested before.
             */
            for (site_t wallSiteLocalIndex = 0; wallSiteLocalIndex < wallSitesCount; wallSiteLocalIndex++)
            {
              site_t streamedToSite = firstWallSite + wallSiteLocalIndex;
              const geometry::Site<geometry::LatticeData> streamedSite =
                  latDat->GetSite(streamedToSite);
              distribn_t* streamedToFNew = latDat->GetFNew(lb::lattices::D3Q15::NUMVECTORS
                  * streamedToSite);

              for (unsigned int streamedDirection = 0; streamedDirection
                  < lb::lattices::D3Q15::NUMVECTORS; ++streamedDirection)
              {
                unsigned oppDirection = lb::lattices::D3Q15::INVERSEDIRECTIONS[streamedDirection];

                // Index of the site streaming to streamedToSite via direction streamedDirection
                site_t streamerIndex =
                    streamedSite.GetStreamedIndex<lb::lattices::D3Q15> (oppDirection);

                // Is streamerIndex a valid index?
                if (streamerIndex >= 0 && streamerIndex < (lb::lattices::D3Q15::NUMVECTORS
                    * latDat->GetLocalFluidSiteCount()))
                {
                  // The streamer index is a valid index in the domain, therefore stream and collide has happened
                  site_t streamerSiteId = streamerIndex / lb::lattices::D3Q15::NUMVECTORS;

                  // Calculate streamerFOld at this site.
                  distribn_t streamerFOld[lb::lattices::D3Q15::NUMVECTORS];
                  LbTestsHelper::InitialiseAnisotropicTestData<lb::lattices::D3Q15>(streamerSiteId,
                                                                                    streamerFOld);

                  // Calculate what the value streamed to site streamedToSite should be.
                  lb::kernels::HydroVars<lb::kernels::LBGK<lb::lattices::D3Q15> >
                      streamerHydroVars(streamerFOld);
                  normalCollision->CalculatePreCollision(streamerHydroVars, streamedSite);

                  normalCollision->Collide(lbmParams, streamerHydroVars);

                  // F_new should be equal to the value that was streamed from this other site
                  // in the same direction as we're streaming from.
                  CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("SimpleCollideAndStream, StreamAndCollide",
                                                       streamerHydroVars.GetFPostCollision()[streamedDirection],
                                                       streamedToFNew[streamedDirection],
                                                       allowedError);
                }
                else
                {
                  // The streamer index shows that no one has streamed to streamedToSite direction
                  // streamedDirection, therefore bounce back has happened in that site for that direction

                  // Initialise streamedToSiteFOld with the original data
                  distribn_t streamerToSiteFOld[lb::lattices::D3Q15::NUMVECTORS];
                  LbTestsHelper::InitialiseAnisotropicTestData<lb::lattices::D3Q15>(streamedToSite,
                                                                                    streamerToSiteFOld);
                  lb::kernels::HydroVars<lb::kernels::LBGK<lb::lattices::D3Q15> >
                      hydroVars(streamerToSiteFOld);
                  normalCollision->CalculatePreCollision(hydroVars, streamedSite);

                  // Simulate post-collision using the collision operator.
                  normalCollision->Collide(lbmParams, hydroVars);

                  // After streaming FNew in a given direction must be f post-collision in the opposite direction
                  // following collision
                  std::stringstream msg;
                  msg << "Junk&Yang bounce-back equivalent: site " << streamedToSite
                      << " direction " << streamedDirection;
                  CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(msg.str(),
                                                       streamedToFNew[streamedDirection],
                                                       hydroVars.GetFPostCollision()[oppDirection],
                                                       allowedError);
                }
              }
            }
          }

          void TestNashZerothOrderPressureIolet()
          {
            lb::iolets::BoundaryValues inletBoundary(geometry::INLET_TYPE,
                                                     latDat,
                                                     simConfig->GetInlets(),
                                                     simState,
                                                     Comms(),
                                                     *unitConverter);

            initParams.boundaryObject = &inletBoundary;

            lb::streamers::NashZerothOrderPressureIolet<lb::collisions::Normal<lb::kernels::LBGK<
                lb::lattices::D3Q15> > >::Type ioletCollider(initParams);

            for (double assignedWallDistance = 0.4; assignedWallDistance < 1.0; assignedWallDistance
                += 0.5)
            {
              // Initialise fOld in the lattice data. We choose values so that each site has
              // an anisotropic distribution function, and that each site's function is
              // distinguishable.
              LbTestsHelper::InitialiseAnisotropicTestData<lb::lattices::D3Q15>(latDat);

              // Make some fairly arbitrary choices early on.
              const site_t chosenSite = 0;
              const int chosenBoundaryId = 0;
              const geometry::Site<geometry::LatticeData>& streamer = latDat->GetSite(chosenSite);

              const Direction chosenUnstreamedDirection = 5;
              const Direction chosenIoletDirection =
                  lb::lattices::D3Q15::INVERSEDIRECTIONS[chosenUnstreamedDirection];
              const util::Vector3D<distribn_t> ioletNormal =
                  inletBoundary.GetLocalIolet(chosenBoundaryId)->GetNormal();

              // Enforce that there's a boundary in the iolet direction.
              latDat->SetHasIolet(chosenSite, chosenIoletDirection);
              latDat->SetBoundaryDistance(chosenSite, chosenIoletDirection, assignedWallDistance);
              latDat->SetBoundaryNormal(chosenSite, ioletNormal);
              latDat->SetIoletId(chosenSite, chosenBoundaryId);

              // Perform the collision and streaming.
              ioletCollider.StreamAndCollide<false> (chosenSite,
                                                     1,
                                                     lbmParams,
                                                     latDat,
                                                     *propertyCache);

              // Check each streamed direction.
              for (Direction streamedDirection = 0; streamedDirection
                  < lb::lattices::D3Q15::NUMVECTORS; ++streamedDirection)
              {
                // Calculate the distributions at the chosen site up to post-collision.
                distribn_t streamerFOld[lb::lattices::D3Q15::NUMVECTORS];
                LbTestsHelper::InitialiseAnisotropicTestData<lb::lattices::D3Q15>(chosenSite,
                                                                                  streamerFOld);

                lb::kernels::HydroVars<lb::kernels::LBGK<lb::lattices::D3Q15> >
                    streamerHydroVars(streamerFOld);
                normalCollision->CalculatePreCollision(streamerHydroVars, streamer);
                normalCollision->Collide(lbmParams, streamerHydroVars);

                // Calculate the streamed-to index.
                const site_t streamedIndex =
                    streamer.GetStreamedIndex<lb::lattices::D3Q15> (streamedDirection);

                // Check that simple collide and stream has happened when appropriate.
                // Is streamerIndex a valid index? (And is it not in one of the directions
                // that has been meddled with for the test)?
                if (!streamer.HasIolet(streamedDirection) && streamedIndex >= 0 && streamedIndex
                    < (lb::lattices::D3Q15::NUMVECTORS * latDat->GetLocalFluidSiteCount()))
                {
                  distribn_t streamedToFNew = *latDat->GetFNew(streamedIndex);

                  // F_new should be equal to the value that was streamed from this other site
                  // in the same direction as we're streaming from.
                  CPPUNIT_ASSERT_DOUBLES_EQUAL(streamerHydroVars.GetFPostCollision()[streamedDirection],
                                               streamedToFNew,
                                               allowedError);
                }

                // Next, handle the case where this is the direction where we're checking for
                // behaviour with a wall. I.e. are we correctly filling distributions that aren't
                // streamed-to by simple streaming?
                if (streamedDirection == chosenUnstreamedDirection)
                {
                  // The streamer works by assuming the presence of a 'ghost' site, just beyond the
                  // iolet. The density of the ghost site is extrapolated from the iolet density
                  // and the density of the fluid site.
                  distribn_t ghostSiteDensity = inletBoundary.GetBoundaryDensity(chosenBoundaryId);

                  // The velocity of the ghost site is the component of the fluid site's velocity
                  // along the iolet normal.
                  util::Vector3D<distribn_t> ghostSiteVelocity = ioletNormal
                      * (streamerHydroVars.momentum / streamerHydroVars.density).Dot(ioletNormal);

                  util::Vector3D<distribn_t> ghostSiteMomentum = ghostSiteVelocity
                      * ghostSiteDensity;

                  distribn_t ghostPostCollision[lb::lattices::D3Q15::NUMVECTORS];

                  LbTestsHelper::CalculateLBGKEqmF<lb::lattices::D3Q15>(ghostSiteDensity,
                                                                        ghostSiteMomentum.x,
                                                                        ghostSiteMomentum.y,
                                                                        ghostSiteMomentum.z,
                                                                        ghostPostCollision);

                  CPPUNIT_ASSERT_DOUBLES_EQUAL(latDat->GetFNew(chosenSite
                                                   * lb::lattices::D3Q15::NUMVECTORS)[chosenUnstreamedDirection],
                                               ghostPostCollision[chosenUnstreamedDirection],
                                               allowedError);
                }
              }
            }
          }

          void TestNashZerothOrderPressureBB()
          {
            lb::iolets::BoundaryValues inletBoundary(geometry::INLET_TYPE,
                                                     latDat,
                                                     simConfig->GetInlets(),
                                                     simState,
                                                     Comms(),
                                                     *unitConverter);

            initParams.boundaryObject = &inletBoundary;

            lb::streamers::NashZerothOrderPressureIoletSBB<lb::collisions::Normal<
                lb::kernels::LBGK<lb::lattices::D3Q15> > >::Type ioletCollider(initParams);

            for (double assignedIoletDistance = 0.4; assignedIoletDistance < 1.0; assignedIoletDistance
                += 0.5)
            {
              // Initialise fOld in the lattice data. We choose values so that each site has
              // an anisotropic distribution function, and that each site's function is
              // distinguishable.
              LbTestsHelper::InitialiseAnisotropicTestData<lb::lattices::D3Q15>(latDat);

              // Make some fairly arbitrary choices early on.
              const site_t chosenSite = 0;
              const int chosenBoundaryId = 0;
              const geometry::Site<geometry::LatticeData>& streamer = latDat->GetSite(chosenSite);

              const Direction chosenWallDirection = 11;
              const Direction chosenUnstreamedDirection = 5;
              const Direction chosenIoletDirection =
                  lb::lattices::D3Q15::INVERSEDIRECTIONS[chosenUnstreamedDirection];
              const util::Vector3D<distribn_t> ioletNormal =
                  inletBoundary.GetLocalIolet(chosenBoundaryId)->GetNormal();

              // Enforce that there's a boundary in the iolet direction.
              latDat->SetHasIolet(chosenSite, chosenIoletDirection);
              latDat->SetHasWall(chosenSite, chosenWallDirection);
              latDat->SetBoundaryDistance(chosenSite, chosenIoletDirection, assignedIoletDistance);
              latDat->SetBoundaryNormal(chosenSite, ioletNormal);
              latDat->SetIoletId(chosenSite, chosenBoundaryId);

              // Perform the collision and streaming.
              ioletCollider.StreamAndCollide<false> (chosenSite,
                                                     1,
                                                     lbmParams,
                                                     latDat,
                                                     *propertyCache);

              // Check each streamed direction.
              for (Direction streamedDirection = 0; streamedDirection
                  < lb::lattices::D3Q15::NUMVECTORS; ++streamedDirection)
              {
                // Calculate the distributions at the chosen site up to post-collision.
                distribn_t streamerFOld[lb::lattices::D3Q15::NUMVECTORS];
                LbTestsHelper::InitialiseAnisotropicTestData<lb::lattices::D3Q15>(chosenSite,
                                                                                  streamerFOld);

                lb::kernels::HydroVars<lb::kernels::LBGK<lb::lattices::D3Q15> >
                    streamerHydroVars(streamerFOld);
                normalCollision->CalculatePreCollision(streamerHydroVars, streamer);
                normalCollision->Collide(lbmParams, streamerHydroVars);

                // Calculate the streamed-to index.
                const site_t streamedIndex =
                    streamer.GetStreamedIndex<lb::lattices::D3Q15> (streamedDirection);

                // Check that simple collide and stream has happened when appropriate.
                // Is streamerIndex a valid index? (And is it not in one of the directions
                // that has been meddled with for the test)?
                if (!streamer.HasIolet(streamedDirection) && !streamer.HasWall(streamedDirection)
                    && streamedIndex >= 0 && streamedIndex < (lb::lattices::D3Q15::NUMVECTORS
                    * latDat->GetLocalFluidSiteCount()))
                {
                  distribn_t streamedToFNew = *latDat->GetFNew(streamedIndex);

                  // F_new should be equal to the value that was streamed from this other site
                  // in the same direction as we're streaming from.
                  CPPUNIT_ASSERT_DOUBLES_EQUAL(streamerHydroVars.GetFPostCollision()[streamedDirection],
                                               streamedToFNew,
                                               allowedError);
                }

                Direction inverseDirection =
                    lb::lattices::D3Q15::INVERSEDIRECTIONS[streamedDirection];

                // Check the case by a wall.
                if (streamer.HasWall(streamedDirection))
                {
                  distribn_t streamedToFNew = *latDat->GetFNew(lb::lattices::D3Q15::NUMVECTORS
                      * chosenSite + inverseDirection);

                  CPPUNIT_ASSERT_DOUBLES_EQUAL(streamerHydroVars.GetFPostCollision()[streamedDirection],
                                               streamedToFNew,
                                               allowedError);
                }

                // Next, handle the case where this is the direction where we're checking for
                // behaviour with a iolet. I.e. are we correctly filling distributions that aren't
                // streamed-to by simple streaming?
                if (streamedDirection == chosenUnstreamedDirection)
                {
                  // The streamer works by assuming the presence of a 'ghost' site, just beyond the
                  // iolet. The density of the ghost site is extrapolated from the iolet density
                  // and the density of the fluid site.
                  distribn_t ghostSiteDensity = inletBoundary.GetBoundaryDensity(chosenBoundaryId);

                  // The velocity of the ghost site is the component of the fluid site's velocity
                  // along the iolet normal.
                  util::Vector3D<distribn_t> ghostSiteVelocity = ioletNormal
                      * (streamerHydroVars.momentum / streamerHydroVars.density).Dot(ioletNormal);

                  util::Vector3D<distribn_t> ghostSiteMomentum = ghostSiteVelocity
                      * ghostSiteDensity;

                  distribn_t ghostPostCollision[lb::lattices::D3Q15::NUMVECTORS];

                  LbTestsHelper::CalculateLBGKEqmF<lb::lattices::D3Q15>(ghostSiteDensity,
                                                                        ghostSiteMomentum.x,
                                                                        ghostSiteMomentum.y,
                                                                        ghostSiteMomentum.z,
                                                                        ghostPostCollision);

                  CPPUNIT_ASSERT_DOUBLES_EQUAL(latDat->GetFNew(chosenSite
                                                   * lb::lattices::D3Q15::NUMVECTORS)[chosenUnstreamedDirection],
                                               ghostPostCollision[chosenUnstreamedDirection],
                                               allowedError);
                }
              }
            }
          }

        private:
          lb::MacroscopicPropertyCache* propertyCache;
          lb::collisions::Normal<lb::kernels::LBGK<lb::lattices::D3Q15> >* normalCollision;
      };
      CPPUNIT_TEST_SUITE_REGISTRATION ( StreamerTests);
    }
  }
}

#endif /* HEMELB_UNITTESTS_LBTESTS_STREAMERTESTS_H */
