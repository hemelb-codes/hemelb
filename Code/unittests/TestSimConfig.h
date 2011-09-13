#ifndef HEMELB_UNITTESTS_TESTSIMCONFIG_H
#define HEMELB_UNITTESTS_TESTSIMCONFIG_H

#include "SimConfig.h"

namespace hemelb
{
  namespace unittests
  {
    class TestSimConfig : public SimConfig
    {
      public:
        TestSimConfig() :
          SimConfig()
        {
          InOutLet inlet;
          inlet.PAmp = 1.0;
          inlet.PMax = 81.0;
          inlet.PMin = 79.0;
          inlet.PMean = 80.0;
          inlet.PPhase = 180.0;

          Inlets.push_back(inlet);

          InOutLet outlet;
          outlet.PAmp = 0.0;
          outlet.PMax = 80.0;
          outlet.PMin = 80.0;
          outlet.PMean = 80.0;
          outlet.PPhase = 0.0;

          Outlets.push_back(outlet);

          NumCycles = 10;
          StepsPerCycle = 1000;
        }
    };
  }
}

#endif /* HEMELB_UNITTESTS_TESTSIMCONFIG_H */
