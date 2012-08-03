// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_STEERING_BASIC_SIMULATIONPARAMETERS_H
#define HEMELB_STEERING_BASIC_SIMULATIONPARAMETERS_H

#include "io/writers/xdr/XdrMemWriter.h"
#include "lb/SimulationState.h"

namespace hemelb
{
  namespace steering
  {

    class SimulationParameters
    {

      public:

        int timeStep;
        double time;
        int nInlets;
        double mousePressure;
        double mouseStress;

        static const u_int paramsSizeB = 3 * sizeof(int) + 3 * sizeof(double);

        SimulationParameters();
        ~SimulationParameters();
        void pack(io::writers::xdr::XdrWriter* writer);

      private:

    };

  }
}

#endif // HEMELB_STEERING_BASIC_SIMULATIONPARAMETERS_H
