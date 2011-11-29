#include "io/writers/xdr/XdrMemWriter.h"

#include "vis/Control.h"
#include "steering/basic/SimulationParameters.h"

namespace hemelb
{
  namespace steering
  {

    SimulationParameters::SimulationParameters()
    {
      // C'tor initialises to the following defaults.

      timeStep = 0;
      time = 0.0;
      cycle = 0;
      nInlets = 3;
      mousePressure = -1.0;
      mouseStress = -1.0;
    }

    SimulationParameters::~SimulationParameters()
    {
    }

    void SimulationParameters::pack(io::writers::xdr::XdrWriter* writer)
    {
      writer->operator <<(timeStep);

      writer->operator <<(time);

      writer->operator <<(cycle);
      writer->operator <<(nInlets);

      writer->operator <<(mousePressure);
      writer->operator <<(mouseStress);
    }

  }
}
