#include "io/XdrMemWriter.h"

#include "vis/Control.h"
#include "steering/basic/SimulationParameters.h"

namespace hemelb
{
  namespace steering
  {

    SimulationParameters::SimulationParameters() :
      paramWriter(params, paramsSizeB)
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

    char* SimulationParameters::pack()
    {
      paramWriter << timeStep;

      paramWriter << time;

      paramWriter << cycle;
      paramWriter << nInlets;

      paramWriter << mousePressure;
      paramWriter << mouseStress;

      return params;
    }

  }
}
