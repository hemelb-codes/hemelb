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
    SteeringComponent::SteeringComponent(Network* network,
                                         vis::Control* iVisControl,
                                         steering::ImageSendComponent* imageSendComponent,
                                         net::Net * iNet,
                                         lb::SimulationState * iSimState,
                                         configuration::SimConfig* iSimConfig,
                                         util::UnitConverter* iUnits) :
        net::PhasedBroadcastRegular<false, 1, 0, true, false>(iNet, iSimState, SPREADFACTOR),
        mSimState(iSimState), mVisControl(iVisControl), imageSendComponent(imageSendComponent), mUnits(iUnits)
    {
      Reset(iSimConfig);
      AssignValues();
    }

    bool SteeringComponent::RequiresSeparateSteeringCore()
    {
      return false;
    }

    void SteeringComponent::ProgressFromParent(unsigned long splayNumber)
    {

    }
    void SteeringComponent::ProgressToChildren(unsigned long splayNumber)
    {

    }

    void SteeringComponent::TopNodeAction()
    {

    }

    void SteeringComponent::Effect()
    {
      AssignValues();
    }
  }
}
