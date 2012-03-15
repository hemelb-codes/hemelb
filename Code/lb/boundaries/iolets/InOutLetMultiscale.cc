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
        void InOutLetMultiscale::DoIO(TiXmlElement *iParent, bool iIsLoading, configuration::SimConfig* iSimConfig)
        {
          iSimConfig->DoIOForMultiscaleInOutlet(iParent, iIsLoading, this);
        }
      }
    }
  }
}
