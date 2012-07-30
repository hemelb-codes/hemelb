#ifndef HEMELB_UNITTESTS_ONEINONEOUTSIMCONFIG_H
#define HEMELB_UNITTESTS_ONEINONEOUTSIMCONFIG_H

#include "configuration/SimConfig.h"

namespace hemelb
{
  namespace unittests
  {
    class OneInOneOutSimConfig : public configuration::SimConfig
    {
      public:
        OneInOneOutSimConfig() :
          configuration::SimConfig()
        {
          lb::boundaries::iolets::InOutLetCosine* inlet = new lb::boundaries::iolets::InOutLetCosine();
          inlet->SetPressureAmp(1.0);
          inlet->SetPressureMean(80.0);
          inlet->SetPhase(PI);
          inlet->SetPeriod(60.0/70.0);
          inlet->SetNormal(util::Vector3D<float>(-3,4,-9));

          inlets.push_back(inlet);

          lb::boundaries::iolets::InOutLetCosine* outlet = new lb::boundaries::iolets::InOutLetCosine();
          outlet->SetPressureAmp(0.0);
          outlet->SetPressureMean(80.0);
          outlet->SetPhase(0.0);
          outlet->SetNormal(util::Vector3D<float>(2,-1,4));

          outlets.push_back(outlet);

          totalTimeSteps = 10000;
          timeStepLength = 60.0/(70.0*1000.0);
        }
    };
  }
}

#endif /* HEMELB_UNITTESTS_ONEINONEOUTSIMCONFIG_H */
