#ifndef HEMELB_MULTISCALE_MULTISCALE_SIMULATION_MASTER_H
#define HEMELB_MULTISCALE_MULTISCALE_SIMULATION_MASTER_H
#include <vector>
#include "multiscale/Intercommunicator.h"
#include "SimulationMaster.h"

namespace hemelb
{
  namespace multiscale
  {
    template<class Intercommunicator> class MultiscaleSimulationMaster : public SimulationMaster
    {
      public:
        MultiscaleSimulationMaster(hemelb::configuration::CommandLine &options, Intercommunicator & aintercomms) :
            SimulationMaster(options), intercomms(aintercomms), multiscaleIoletType("inoutlet")
        {
          lb::boundaries::iolets::InOutLetMultiscale::DefineType(multiscaleIoletType);
          for (unsigned int i = 0; i < inletValues->GetLocalIoletCount(); i++)
          {
            if (inletValues->GetLocalIolet(i)->IsRegistrationRequired())
            {
              static_cast<lb::boundaries::iolets::InOutLetMultiscale*>(inletValues->GetLocalIolet(i))->Register(intercomms,multiscaleIoletType);
            }
          }
        }

        void DoTimeStep()
        {
          if (!intercomms.ShouldAdvance())
          {
            hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("HemeLB waiting pending multiscale siblings.");
            return;
          }
          intercomms.GetFromMultiscale();
          SimulationMaster::DoTimeStep();
          hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("HemeLB advanced to time %f.",
                                                                              GetState()->GetTime());
          intercomms.AdvanceTime(GetState()->GetTime());
          intercomms.SendToMultiscale();
        }
        Intercommunicator &intercomms;
        typename Intercommunicator::IntercommunicandTypeT multiscaleIoletType;
    };
  }
}

#endif // HEMELB_MULTISCALE_MULTISCALE_SIMULATION_MASTER_H
