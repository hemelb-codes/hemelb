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
      class StreamerTests : public CppUnit::TestFixture
      {
        public:
          void setUp()
          {
            int args = 1;
            char** argv = NULL;
            bool success;
            topology::NetworkTopology::Instance()->Init(&args, &argv, &success);

            latDat = new TestLatticeData();
            simConfig = new TestSimConfig();
            simState = new lb::SimulationState(simConfig->StepsPerCycle, simConfig->NumCycles);
            lbmParams = new lb::LbmParameters(PULSATILE_PERIOD_s
                                                  / (distribn_t) simState->GetTimeStepsPerCycle(),
                                              latDat->GetVoxelSize());
            unitConverter = new util::UnitConverter(lbmParams, simState, latDat);

            // Initialise the collision.
            lb::kernels::InitParams initParams;
            initParams.latDat = latDat;
            normalCollision = new lb::collisions::Normal<lb::kernels::LBGK>(initParams);
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
              distribn_t* fOld = latDat->GetFOld(ii);
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
            for (site_t ii = 0; ii < latDat->GetLocalFluidSiteCount(); ++ii)
            {
              distribn_t* fNew = latDat->GetFOld(ii);
              for (unsigned int jj = 0; jj < D3Q15::NUMVECTORS; ++jj)
              {
                site_t streamedIndex = latDat->GetStreamedIndex(ii, jj);

                // If this site streamed somewhere sensible, it must have been streamed to.
                if (streamedIndex >= 0
                    && streamedIndex < (D3Q15::NUMVECTORS * latDat->GetLocalFluidSiteCount()))
                {
                  site_t streamedSite = streamedIndex / D3Q15::NUMVECTORS;

                  // Calculate fOld at this site.
                  distribn_t fOld[D3Q15::NUMVECTORS];
                  for (unsigned int kk = 0; kk < D3Q15::NUMVECTORS; ++kk)
                  {
                    fOld[kk] = ((distribn_t) (kk + 1)) / 10.0
                        + ((distribn_t) (streamedSite + 1)) / 100.0;
                  }

                  // Calculate what the value streamed to site ii should be.
                  lb::kernels::HydroVars<lb::kernels::LBGK> hydroVars(fOld);
                  normalCollision->CalculatePreCollision(hydroVars, streamedSite);
                  for (unsigned int kk = 0; kk < D3Q15::NUMVECTORS; ++kk)
                  {
                    hydroVars.f_neq[kk] = hydroVars.f[kk] - hydroVars.f_eq[kk];
                  }

                  // F_new should be equal to the value that was streamed from this other site
                  // in the direction /opposite/ to the direction we're streaming from.
                  CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("SimpleCollideAndStream, StreamAndCollide",
                                                       normalCollision->Collide(lbmParams,
                                                                                D3Q15::INVERSEDIRECTIONS[jj],
                                                                                hydroVars),
                                                       fNew[jj],
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
