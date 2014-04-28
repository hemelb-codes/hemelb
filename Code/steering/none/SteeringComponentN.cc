
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "steering/SteeringComponent.h"

namespace hemelb
{
  namespace steering
  {

    /**
     * In the 'no steering' version of the steering controller, we don't
     * need to override any of the base classes methods, because we don't
     * need to do anything.
     *
     * @param iNet
     * @param iNetTop
     * @param iSimState
     * @return
     */
    SteeringComponent::SteeringComponent(Network* network, vis::Control* iVisControl,
                                         steering::ImageSendComponent* imageSendComponent,
                                         net::Net * iNet, lb::SimulationState * iSimState,
                                         configuration::SimConfig* iSimConfig,
                                         const util::UnitConverter* iUnits,
                                         reporting::Timers& timings) :
        net::CollectiveAction(iNet->GetCommunicator(), timings), mSimState(iSimState),
            mVisControl(iVisControl), imageSendComponent(imageSendComponent), mUnits(iUnits),
            simConfig(iSimConfig), privateSteeringParams(STEERABLE_PARAMETERS + 1)
    {
      ClearValues();
      AssignValues();
    }

    bool SteeringComponent::RequiresSeparateSteeringCore()
    {
      return false;
    }

    void SteeringComponent::PreSend()
    {
    }
    void SteeringComponent::Send()
    {
    }
    void SteeringComponent::PostReceive()
    {
      AssignValues();
    }
  }
}
