// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

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

          inlets.push_back(inlet);

          lb::boundaries::iolets::InOutLetCosine* outlet = new lb::boundaries::iolets::InOutLetCosine();
          outlet->SetPressureAmp(0.0);
          outlet->SetPressureMean(80.0);
          outlet->SetPhase(0.0);

          outlets.push_back(outlet);

          totalTimeSteps = 10000;
          timeStepLength = 60.0/(70.0*1000.0);
        }
    };
  }
}

#endif /* HEMELB_UNITTESTS_ONEINONEOUTSIMCONFIG_H */
