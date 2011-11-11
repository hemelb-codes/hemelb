#ifndef HEMELB_UNITTESTS_LBTESTS_STREAMERTESTS_H
#define HEMELB_UNITTESTS_LBTESTS_STREAMERTESTS_H

#include <cppunit/TestFixture.h>
#include <iostream>
#include <sstream>

#include "lb/streamers/Streamers.h"

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
      class StreamerTests : public CppUnit::TestFixture
      {
        public:

          void setUp()
          {
            int args = 1;
            char** argv = NULL;
            bool success;
            topology::NetworkTopology::Instance()->Init(&args, &argv, &success);

            latDat = new FourCubeLatticeData();
            simConfig = new OneInOneOutSimConfig();
            simState = new lb::SimulationState(simConfig->StepsPerCycle, simConfig->NumCycles);
            lbmParams = new lb::LbmParameters(PULSATILE_PERIOD_s
                                                  / (distribn_t) simState->GetTimeStepsPerCycle(),
                                              latDat->GetVoxelSize());
            unitConverter = new util::UnitConverter(lbmParams, simState, latDat->GetVoxelSize());

            // Initialise the collision.
            lb::kernels::InitParams initParams;
            initParams.latDat = latDat;
            normalCollision = new lb::collisions::Normal<lb::kernels::LBGK>(initParams);

            simpleCollideAndStream = new lb::streamers::SimpleCollideAndStream<
                lb::collisions::Normal<lb::kernels::LBGK> >(initParams);
          }

          void tearDown()
          {
            delete latDat;
            delete simConfig;
            delete simState;
            delete unitConverter;
            delete lbmParams;

            delete normalCollision;

            delete simpleCollideAndStream;
            delete fInterpolation;
            delete regularised;
            delete guoZhengShi;
            delete simpleBounceBack;
          }

          void TestSimpleCollideAndStream()
          {
            // Initialise fOld in the lattice data. We choose values so that each site has
            // an anisotropic distribution function, and that each site's function is
            // distinguishable.
            LbTestsHelper::InitialiseAnisotropicTestData(latDat);

            // Use the streaming operator on the entire lattice.
            simpleCollideAndStream->StreamAndCollide<false> (0,
                                                             latDat->GetLocalFluidSiteCount(),
                                                             lbmParams,
                                                             latDat,
                                                             NULL);

            // Now, go over each lattice site and check each value in f_new is correct.
            for (site_t streamedToSite = 0; streamedToSite < latDat->GetLocalFluidSiteCount(); ++streamedToSite)
            {
              distribn_t* streamedToFNew = latDat->GetFNew(D3Q15::NUMVECTORS * streamedToSite);

              for (unsigned int streamedDirection = 0; streamedDirection < D3Q15::NUMVECTORS; ++streamedDirection)
              {

                site_t streamerIndex =
                    latDat->GetStreamedIndex(streamedToSite,
                                             D3Q15::INVERSEDIRECTIONS[streamedDirection]);

                // If this site streamed somewhere sensible, it must have been streamed to.
                if (streamerIndex >= 0 && streamerIndex < (D3Q15::NUMVECTORS
                    * latDat->GetLocalFluidSiteCount()))
                {
                  site_t streamerSiteId = streamerIndex / D3Q15::NUMVECTORS;

                  // Calculate streamerFOld at this site.
                  distribn_t streamerFOld[D3Q15::NUMVECTORS];
                  LbTestsHelper::InitialiseAnisotropicTestData(streamerSiteId, streamerFOld);

                  // Calculate what the value streamed to site streamedToSite should be.
                  lb::kernels::HydroVars<lb::kernels::LBGK> streamerHydroVars(streamerFOld);
                  normalCollision->CalculatePreCollision(streamerHydroVars, streamerSiteId);

                  // F_new should be equal to the value that was streamed from this other site
                  // in the same direction as we're streaming from.
                  CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("SimpleCollideAndStream, StreamAndCollide",
                                                       normalCollision->Collide(lbmParams,
                                                                                streamedDirection,
                                                                                streamerHydroVars),
                                                       streamedToFNew[streamedDirection],
                                                       allowedError);
                }
              }
            }
          }

          void TestFInterpolation()
          {
            // Initialise fOld in the lattice data. We choose values so that each site has
            // an anisotropic distribution function, and that each site's function is
            // distinguishable.
            LbTestsHelper::InitialiseAnisotropicTestData(latDat);

            fInterpolation->StreamAndCollide<false> (0,
                                                     latDat->GetLocalFluidSiteCount(),
                                                     lbmParams,
                                                     latDat,
                                                     NULL);
            fInterpolation->PostStep<false> (0,
                                             latDat->GetLocalFluidSiteCount(),
                                             lbmParams,
                                             latDat,
                                             NULL);

            // Now, go over each lattice site and check each value in f_new is correct.
            for (site_t streamedToSite = 0; streamedToSite < latDat->GetLocalFluidSiteCount(); ++streamedToSite)
            {
              distribn_t* streamedToFNew = latDat->GetFNew(D3Q15::NUMVECTORS * streamedToSite);

              for (unsigned int streamedDirection = 0; streamedDirection < D3Q15::NUMVECTORS; ++streamedDirection)
              {
                unsigned int oppDirection = D3Q15::INVERSEDIRECTIONS[streamedDirection];

                site_t streamerIndex = latDat->GetStreamedIndex(streamedToSite, oppDirection);

                // If this site streamed somewhere sensible, it must have been streamed to.
                if (streamerIndex >= 0 && streamerIndex < (D3Q15::NUMVECTORS
                    * latDat->GetLocalFluidSiteCount()))
                {
                  site_t streamerSiteId = streamerIndex / D3Q15::NUMVECTORS;

                  // Calculate streamerFOld at this site.
                  distribn_t streamerFOld[D3Q15::NUMVECTORS];
                  LbTestsHelper::InitialiseAnisotropicTestData(streamerSiteId, streamerFOld);

                  // Calculate what the value streamed to site streamedToSite should be.
                  lb::kernels::HydroVars<lb::kernels::LBGK> streamerHydroVars(streamerFOld);
                  normalCollision->CalculatePreCollision(streamerHydroVars, streamerSiteId);

                  // F_new should be equal to the value that was streamed from this other site
                  // in the same direction as we're streaming from.
                  CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("SimpleCollideAndStream, StreamAndCollide",
                                                       normalCollision->Collide(lbmParams,
                                                                                streamedDirection,
                                                                                streamerHydroVars),
                                                       streamedToFNew[streamedDirection],
                                                       allowedError);
                }

                else
                {
                  CPPUNIT_ASSERT_MESSAGE("Expected to find a boundary"
                                           "opposite an unstreamed-to direction",
                                         latDat->HasBoundary(streamedToSite, oppDirection));

                  // To verify the operation of the f-interpolation boundary condition, we'll need:
                  // - the distance to the wall * 2
                  distribn_t twoQ = 2.0 * latDat->GetCutDistance(streamedToSite, oppDirection);

                  // - the post-collision distribution at the current site.
                  distribn_t streamedToSitePostColl[D3Q15::NUMVECTORS];

                  // (initialise it to f_old).
                  LbTestsHelper::InitialiseAnisotropicTestData(streamedToSite,
                                                               streamedToSitePostColl);

                  lb::kernels::HydroVars<lb::kernels::LBGK> hydroVars(streamedToSitePostColl);

                  normalCollision->CalculatePreCollision(hydroVars, streamedToSite);

                  // (find post-collision values using the collision operator).
                  for (unsigned int kk = 0; kk < D3Q15::NUMVECTORS; ++kk)
                  {
                    streamedToSitePostColl[kk] = normalCollision->Collide(lbmParams, kk, hydroVars);
                  }

                  // - finally, the post-collision distribution at the site which is one further
                  // away from the wall in this direction.
                  distribn_t awayFromWallPostColl[D3Q15::NUMVECTORS];

                  site_t awayFromWallIndex = latDat->GetStreamedIndex(streamedToSite,
                                                                      streamedDirection)
                      / D3Q15::NUMVECTORS;

                  // If there's a valid index in that direction, use f-interpolation
                  if (awayFromWallIndex >= 0 && awayFromWallIndex
                      < latDat->GetLocalFluidSiteCount())
                  {

                    // (initialise it to f_old).
                    LbTestsHelper::InitialiseAnisotropicTestData(awayFromWallIndex,
                                                                 awayFromWallPostColl);

                    lb::kernels::HydroVars<lb::kernels::LBGK>
                        awayFromWallsHydroVars(awayFromWallPostColl);

                    normalCollision->CalculatePreCollision(awayFromWallsHydroVars,
                                                           awayFromWallIndex);

                    // (find post-collision values using the collision operator).
                    for (unsigned int kk = 0; kk < D3Q15::NUMVECTORS; ++kk)
                    {
                      awayFromWallPostColl[kk] = normalCollision->Collide(lbmParams,
                                                                          kk,
                                                                          awayFromWallsHydroVars);
                    }

                    distribn_t toWallOld = streamedToSitePostColl[oppDirection];

                    distribn_t toWallNew = awayFromWallPostColl[oppDirection];

                    distribn_t oppWallOld = streamedToSitePostColl[streamedDirection];

                    // The streamed value should be as given below.
                    distribn_t streamed = (twoQ < 1.0)
                      ? (toWallNew + twoQ * (toWallOld - toWallNew))
                      : (oppWallOld + (1. / twoQ) * (toWallOld - oppWallOld));

                    std::stringstream msg(std::stringstream::in);
                    msg << "FInterpolation, PostStep: site " << streamedToSite << " direction "
                        << streamedDirection;

                    // Assert that this is the case.
                    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(msg.str(),
                                                         streamed,
                                                         * (latDat->GetFNew(streamedToSite
                                                             * D3Q15::NUMVECTORS
                                                             + streamedDirection)),
                                                         allowedError);
                  }

                  // With no valid lattice site, simple bounce-back will be performed.
                  else
                  {
                    std::stringstream msg(std::stringstream::in);
                    msg << "FInterpolation, PostStep by simple bounce-back: site "
                        << streamedToSite << " direction " << streamedDirection;

                    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(msg.str(),
                                                         streamedToSitePostColl[oppDirection],
                                                         * (latDat->GetFNew(streamedToSite
                                                             * D3Q15::NUMVECTORS
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
            LbTestsHelper::InitialiseAnisotropicTestData(latDat);

            site_t firstWallSite = latDat->GetInnerCollisionCount(0);
            site_t wallSitesCount = latDat->GetInnerCollisionCount(1) - firstWallSite;

            // Check that the lattice has sites labeled as wall (otherwise this test is void)
            CPPUNIT_ASSERT_EQUAL_MESSAGE("Number of non-wall sites", wallSitesCount, (site_t) 16);

            site_t offset = 0;

            // Mid-Fluid sites use simple collide and stream
            simpleCollideAndStream->StreamAndCollide<false> (offset,
                                                             latDat->GetInnerCollisionCount(0),
                                                             lbmParams,
                                                             latDat,
                                                             NULL);
            offset += latDat->GetInnerCollisionCount(0);

            // Wall sites use simple bounce back
            simpleBounceBack->StreamAndCollide<false> (offset,
                                                       latDat->GetInnerCollisionCount(1),
                                                       lbmParams,
                                                       latDat,
                                                       NULL);
            offset += latDat->GetInnerCollisionCount(1);

            // Consider inlet/outlets and their walls as mid-fluid sites
            simpleCollideAndStream->StreamAndCollide<false> (offset,
                                                             latDat->GetLocalFluidSiteCount()
                                                                 - offset,
                                                             lbmParams,
                                                             latDat,
                                                             NULL);
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
              distribn_t* streamedToFNew = latDat->GetFNew(D3Q15::NUMVECTORS * streamedToSite);

              for (unsigned int streamedDirection = 0; streamedDirection < D3Q15::NUMVECTORS; ++streamedDirection)
              {
                unsigned oppDirection = D3Q15::INVERSEDIRECTIONS[streamedDirection];

                // Index of the site streaming to streamedToSite via direction streamedDirection
                site_t streamerIndex = latDat->GetStreamedIndex(streamedToSite, oppDirection);

                // Is streamerIndex a valid index?
                if (streamerIndex >= 0 && streamerIndex < (D3Q15::NUMVECTORS
                    * latDat->GetLocalFluidSiteCount()))
                {
                  // The streamer index is a valid index in the domain, therefore stream and collide has happened
                  site_t streamerSiteId = streamerIndex / D3Q15::NUMVECTORS;

                  // Calculate streamerFOld at this site.
                  distribn_t streamerFOld[D3Q15::NUMVECTORS];
                  LbTestsHelper::InitialiseAnisotropicTestData(streamerSiteId, streamerFOld);

                  // Calculate what the value streamed to site streamedToSite should be.
                  lb::kernels::HydroVars<lb::kernels::LBGK> streamerHydroVars(streamerFOld);
                  normalCollision->CalculatePreCollision(streamerHydroVars, streamerSiteId);

                  // F_new should be equal to the value that was streamed from this other site
                  // in the same direction as we're streaming from.
                  CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("SimpleCollideAndStream, StreamAndCollide",
                                                       normalCollision->Collide(lbmParams,
                                                                                streamedDirection,
                                                                                streamerHydroVars),
                                                       streamedToFNew[streamedDirection],
                                                       allowedError);
                }
                else
                {
                  // The streamer index shows that no one has streamed to streamedToSite direction
                  // streamedDirection, therefore bounce back has happened in that site for that direction

                  // Initialise streamedToSiteFOld with the original data
                  distribn_t streamerToSiteFOld[D3Q15::NUMVECTORS];
                  LbTestsHelper::InitialiseAnisotropicTestData(streamedToSite, streamerToSiteFOld);
                  lb::kernels::HydroVars<lb::kernels::LBGK> hydroVars(streamerToSiteFOld);
                  normalCollision->CalculatePreCollision(hydroVars, streamedToSite);

                  // Simulate post-collision using the collision operator.
                  distribn_t streamedToSitePostColl[D3Q15::NUMVECTORS];
                  for (unsigned int kk = 0; kk < D3Q15::NUMVECTORS; ++kk)
                  {
                    streamedToSitePostColl[kk] = normalCollision->Collide(lbmParams, kk, hydroVars);
                  }

                  // After streaming FNew in a given direction must be f post-collision in the opposite direction
                  // following collision
                  std::stringstream msg(std::stringstream::in);
                  msg << "Simple bounce-back: site " << streamedToSite << " direction "
                      << streamedDirection;
                  CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(msg.str(),
                                                       streamedToFNew[streamedDirection],
                                                       streamedToSitePostColl[oppDirection],
                                                       allowedError);
                }
              }
            }
          }

        private:
          geometry::LatticeData* latDat;
          SimConfig* simConfig;
          lb::SimulationState* simState;
          util::UnitConverter* unitConverter;
          lb::LbmParameters* lbmParams;

          lb::collisions::Normal<lb::kernels::LBGK>* normalCollision;

          lb::streamers::SimpleCollideAndStream<lb::collisions::Normal<lb::kernels::LBGK> >
              * simpleCollideAndStream;
          lb::streamers::FInterpolation<lb::collisions::Normal<lb::kernels::LBGK> >
              * fInterpolation;
          lb::streamers::Regularised<lb::collisions::Normal<lb::kernels::LBGK> > * regularised;
          lb::streamers::GuoZhengShi<lb::collisions::Normal<lb::kernels::LBGK> > * guoZhengShi;
          lb::streamers::SimpleBounceBack<lb::collisions::Normal<lb::kernels::LBGK> >
              * simpleBounceBack;

          //          distribn_t allowedError;

      };
    }
  }
}

#endif /* HEMELB_UNITTESTS_LBTESTS_STREAMERTESTS_H */
