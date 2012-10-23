// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_MULTISCALE_MULTISCALESIMULATIONMASTER_H
#define HEMELB_MULTISCALE_MULTISCALESIMULATIONMASTER_H
#include <vector>
#include "multiscale/Intercommunicator.h"
#include "SimulationMaster.h"

namespace hemelb
{
  namespace multiscale
  {
    /***
     * Instead of adding multiscale functionality to the standard simulation master, we keep this here,
     * so the main code can be read without thinking about multiscale.
     */
    template<class Intercommunicator> class MultiscaleSimulationMaster : public SimulationMaster
    {
      public:
        MultiscaleSimulationMaster(hemelb::configuration::CommandLine &options, Intercommunicator & aintercomms) :
            SimulationMaster(options), intercomms(aintercomms), multiscaleIoletType("inoutlet")
        {
          // We only have one shared object type so far, an iolet.
          lb::boundaries::iolets::InOutLetMultiscale::DefineType(multiscaleIoletType);

          // we only want to register those iolets which are needed on this process.
          // Fortunately, the BoundaryValues instance has worked this out for us.
          for (unsigned int i = 0; i < inletValues->GetLocalIoletCount(); i++)
          {
            // could be a if dynamic_cast<> rather than using a castable? virtual method pattern, if we prefer.
            if (inletValues->GetLocalIolet(i)->IsRegistrationRequired())
            {
              static_cast<lb::boundaries::iolets::InOutLetMultiscale*>(inletValues->GetLocalIolet(i))->Register(intercomms,
                                                                                                                multiscaleIoletType);
            }
          }
          for (unsigned int i = 0; i < outletValues->GetLocalIoletCount(); i++)
          {
            if (outletValues->GetLocalIolet(i)->IsRegistrationRequired())
            {
              static_cast<lb::boundaries::iolets::InOutLetMultiscale*>(outletValues->GetLocalIolet(i))->Register(intercomms,
                                                                                                                 multiscaleIoletType);
            }
          }
          intercomms.ShareInitialConditions();
        }

        void DoTimeStep()
        {

          if (intercomms.DoMultiscale(GetState()->GetTime()))
          {
            SimulationMaster::DoTimeStep();
            hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::Singleton>("HemeLB advanced to time %f.",
                                                                           GetState()->GetTime());
          }
          else
          {
            hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::Singleton>("HemeLB waiting pending multiscale siblings.");
            return;
          };
        }
        Intercommunicator &intercomms;
        typename Intercommunicator::IntercommunicandTypeT multiscaleIoletType;
    };
  }
}

#endif // HEMELB_MULTISCALE_MULTISCALE_SIMULATION_MASTER_H
