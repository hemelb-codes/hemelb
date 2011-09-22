#ifndef HEMELB_UNITTESTS_ONEINONEOUTSIMCONFIG_H
#define HEMELB_UNITTESTS_ONEINONEOUTSIMCONFIG_H

#include "SimConfig.h"

namespace hemelb
{
  namespace unittests
  {
    class OneInOneOutSimConfig : public SimConfig
    {
      public:
        OneInOneOutSimConfig() :
          SimConfig()
        {
          lb::boundaries::iolets::InOutLetCosine* inlet = new lb::boundaries::iolets::InOutLetCosine();
          inlet->PressureAmpPhysical = 1.0;
          inlet->PressureMaxPhysical = 81.0;
          inlet->PressureMinPhysical = 79.0;
          inlet->PressureMeanPhysical = 80.0;
          inlet->Phase = 180.0;

          Inlets.push_back(inlet);

          lb::boundaries::iolets::InOutLetCosine* outlet = new lb::boundaries::iolets::InOutLetCosine();
          outlet->PressureAmpPhysical = 0.0;
          outlet->PressureMaxPhysical = 80.0;
          outlet->PressureMinPhysical = 80.0;
          outlet->PressureMeanPhysical = 80.0;
          outlet->Phase = 0.0;

          Outlets.push_back(outlet);

          NumCycles = 10;
          StepsPerCycle = 1000;
        }
    };
  }
}

#endif /* HEMELB_UNITTESTS_ONEINONEOUTSIMCONFIG_H */
