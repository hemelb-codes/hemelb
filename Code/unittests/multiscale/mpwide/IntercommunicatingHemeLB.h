// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_UNITTESTS_MULTISCALE_MPWIDE_INTERCOMMUNICATINGHEMELB_H
#define HEMELB_UNITTESTS_MULTISCALE_MPWIDE_INTERCOMMUNICATINGHEMELB_H

/*
 * IntercommunicatingHemeLB.h
 *
 *  Created on: 26 Jun 2012
 *      Author: derek
 */

#include <map>

namespace hemelb
{
  namespace unittests
  {
    namespace multiscale
    {
      namespace mpwide
      {
        /***
         * Mock intercommunicating entity which looks a bit like a HemeLB conceptually
         * It has an input, and an output, and the flow rate depends on the difference in pressures.
         */
        template<class IntercommuniatorImplementation> class InterCommunicatingHemeLB
        {
          public:
            InterCommunicatingHemeLB(double spaceResolution,
                                     double timeResolution,
                                     std::map<std::string, double> & buffer,
                                     std::map<std::string, bool> &orchestration) :
                inlet(81.0, 0.1), outlet(79.0, 0.1), inOutLetType("inoutlet"), intercomms(buffer, orchestration), timeResolution(timeResolution), spaceResolution(spaceResolution), currentTime(0)
            {
              // The intercommunicators have a shared buffer which represents imaginary communication
              inOutLetType.template RegisterSharedValue<double>("pressure");
              inOutLetType.template RegisterSharedValue<double>("velocity");
              intercomms.RegisterIntercommunicand(inOutLetType, inlet, "boundary1");
              intercomms.RegisterIntercommunicand(inOutLetType, outlet, "boundary2");
              intercomms.ShareInitialConditions();
            }

            void DoLB()
            {
              double resistance = 10.0;
              double velocity = (inlet.GetPressure() - outlet.GetPressure()) / resistance;
              inlet.SetVelocity(velocity);
              outlet.SetVelocity(velocity);
            }

            void Simulate()
            {
              if (intercomms.DoMultiscale(currentTime))
              {
                std::cerr << "MH: " << currentTime << std::endl;
                DoLB();
                currentTime += timeResolution;
              }
            }

            MockIntercommunicand inlet;
            MockIntercommunicand outlet;
            typename IntercommuniatorImplementation::IntercommunicandTypeT inOutLetType;
            IntercommuniatorImplementation intercomms;
            double timeResolution;
            double spaceResolution;
            double currentTime;
        };
      }
    }
  }
}

#endif /* HEMELB_UNITTESTS_MULTISCALE_MPWIDE_INTERCOMMUNICATINGHEMELB_H */
