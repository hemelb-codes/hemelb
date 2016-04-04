
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "steering/SteeringComponent.h"

namespace hemelb
{
  namespace steering
  {
    void SteeringComponent::AssignValues()
    {
      mSimState->SetIsTerminating(1 == (int) privateSteeringParams[SetIsTerminal]);
    }

    void SteeringComponent::ClearValues()
    {
      isConnected = false;

      // signal useful to terminate the simulation
      privateSteeringParams[SetIsTerminal] = 0.0F;

      // Vis_mode
      privateSteeringParams[Mode] = 0.0F;
    }
  }
}
