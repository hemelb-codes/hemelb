#include "lb/boundaries/iolets/InOutLetMultiscale.h"
#include "configuration/SimConfig.h"
#include "topology/NetworkTopology.h"

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
      }
    }
  }
}
