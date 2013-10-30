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
    /**
     * TODO: Figure out what this is supposed to be.
     */
    class OneInOneOutSimConfig : public configuration::SimConfig
    {
      public:
        OneInOneOutSimConfig()
        {
          lb::iolets::InOutLetCosine* inlet = new lb::iolets::InOutLetCosine();
          inlet->SetPressureAmp(1.0);
          inlet->SetPressureMean(80.0);
          inlet->SetPhase(PI);
          inlet->SetPeriod(60.0/70.0);
          inlet->SetNormal(util::Vector3D<Dimensionless>(-3,4,-9));

          inlets.push_back(inlet);

          lb::iolets::InOutLetCosine* outlet = new lb::iolets::InOutLetCosine();
          outlet->SetPressureAmp(0.0);
          outlet->SetPressureMean(80.0);
          outlet->SetPhase(0.0);
          outlet->SetNormal(util::Vector3D<Dimensionless>(2,-1,4));

          outlets.push_back(outlet);

          totalTimeSteps = 10000;
          timeStepLength = 60.0/(70.0*1000.0);
          voxelSizeMetres = 0.01;
          geometryOriginMetres = util::Vector3D<PhysicalDistance>::Zero();
        }
    };
  }
}

#endif /* HEMELB_UNITTESTS_ONEINONEOUTSIMCONFIG_H */
