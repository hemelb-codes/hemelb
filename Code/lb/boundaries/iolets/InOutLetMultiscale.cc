// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include "lb/boundaries/iolets/InOutLetMultiscale.h"
#include "configuration/SimConfig.h"
#include "topology/NetworkTopology.h"
#include "lb/boundaries/BoundaryComms.h"

namespace hemelb
{
  namespace lb
  {
    namespace boundaries
    {
      namespace iolets
      {
        void InOutLetMultiscale::DoIO(TiXmlElement *parent, bool isLoading, configuration::SimConfig* simConfig)
        {
          simConfig->DoIOForMultiscaleInOutlet(parent, isLoading, this);
        }

        /* Distribution of internal pressure values */
        void InOutLetMultiscale::DoComms(bool is_io_proc, LatticeTime time_step)
        {
          hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("DoComms in IoletMultiscale triggered: %s", is_io_proc ? "true" : "false");
          double pressure_array[3];
          //TODO: Change these operators on SharedValue.
          pressure_array[0] = pressure.GetPayload();
          pressure_array[1] = minPressure.GetPayload();
          pressure_array[2] = maxPressure.GetPayload();

          if (is_io_proc)
          {
            //hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::OnePerCore>("comms->SendDoubles().");
            comms->SendDoubles(pressure_array,3);
            //hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::OnePerCore>("comms->FinishSend().");
            }
          else
          {
            //hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::OnePerCore>("comms->ReceiveDoubles().");
            comms->ReceiveDoubles(pressure_array,3);
            pressure.SetPayload(static_cast<PhysicalPressure>(pressure_array[0]));
            minPressure.SetPayload(static_cast<PhysicalPressure>(pressure_array[1]));
            maxPressure.SetPayload(static_cast<PhysicalPressure>(pressure_array[2]));
            hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("Received: %f %f %f", pressure.GetPayload(), minPressure.GetPayload(), maxPressure.GetPayload());
          }
        }
      }
    }
  }
}
