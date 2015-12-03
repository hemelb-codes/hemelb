// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

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
