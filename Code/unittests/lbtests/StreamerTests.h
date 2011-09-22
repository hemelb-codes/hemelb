#ifndef HEMELB_UNITTESTS_LBTESTS_STREAMERTESTS_H
#define HEMELB_UNITTESTS_LBTESTS_STREAMERTESTS_H

#include <cppunit/TestFixture.h>

#include "lb/streamers/Streamers.h"

namespace hemelb
{
  namespace unittests
  {
    namespace lbtests
    {
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
            unitConverter = new util::UnitConverter(lbmParams, simState, latDat);

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
          }

          void TestSimpleCollideAndStream()
          {
            const distribn_t allowedError = 1e-10;

            // Initialise fOld in the lattice data. We choose values so that each site has
            // an anisotropic distribution function, and that each site's function is
            // distinguishable.
            for (site_t ii = 0; ii < latDat->GetLocalFluidSiteCount(); ++ii)
            {
              distribn_t* fOld = latDat->GetFOld(ii * D3Q15::NUMVECTORS);
              for (unsigned int jj = 0; jj < D3Q15::NUMVECTORS; ++jj)
              {
                fOld[jj] = ((distribn_t) (jj + 1)) / 10.0 + ((distribn_t) (ii + 1)) / 100.0;
              }
            }

            // Use the streaming operator on the entire lattice.
            simpleCollideAndStream->StreamAndCollide<false>(0,
                                                            latDat->GetLocalFluidSiteCount(),
                                                            lbmParams,
                                                            latDat,
                                                            NULL);

            // Now, go over each lattice site and check each value in f_new is correct.
            for (site_t streamedToSite = 0; streamedToSite < latDat->GetLocalFluidSiteCount();
                ++streamedToSite)
                {
              distribn_t* streamedToFNew = latDat->GetFNew(D3Q15::NUMVECTORS * streamedToSite);

              for (unsigned int streamedDirection = 0; streamedDirection < D3Q15::NUMVECTORS;
                  ++streamedDirection)
                  {

                site_t streamerIndex =
                    latDat->GetStreamedIndex(streamedToSite,
                                             D3Q15::INVERSEDIRECTIONS[streamedDirection]);

                // If this site streamed somewhere sensible, it must have been streamed to.
                if (streamerIndex >= 0
                    && streamerIndex < (D3Q15::NUMVECTORS * latDat->GetLocalFluidSiteCount()))
                {
                  site_t streamerSiteId = streamerIndex / D3Q15::NUMVECTORS;

                  // Calculate streamerFOld at this site.
                  distribn_t streamerFOld[D3Q15::NUMVECTORS];
                  for (unsigned int kk = 0; kk < D3Q15::NUMVECTORS; ++kk)
                  {
                    streamerFOld[kk] = ((distribn_t) (kk + 1)) / 10.0
                        + ((distribn_t) (streamerSiteId + 1)) / 100.0;
                  }

                  // Calculate what the value streamed to site streamedToSite should be.
                  lb::kernels::HydroVars<lb::kernels::LBGK> streamerHydroVars(streamerFOld);
                  normalCollision->CalculatePreCollision(streamerHydroVars, streamerSiteId);
                  for (unsigned int kk = 0; kk < D3Q15::NUMVECTORS; ++kk)
                  {
                    streamerHydroVars.f_neq[kk] = streamerHydroVars.f[kk]
                        - streamerHydroVars.f_eq[kk];
                  }

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

        private:
          geometry::LatticeData* latDat;
          SimConfig* simConfig;
          lb::SimulationState* simState;
          util::UnitConverter* unitConverter;
          lb::LbmParameters* lbmParams;

          lb::collisions::Normal<lb::kernels::LBGK>* normalCollision;

          lb::streamers::SimpleCollideAndStream<lb::collisions::Normal<lb::kernels::LBGK> >* simpleCollideAndStream;
      };
    }
  }
}

#endif /* HEMELB_UNITTESTS_LBTESTS_STREAMERTESTS_H */
