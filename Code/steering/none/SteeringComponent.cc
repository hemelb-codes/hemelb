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
    SteeringComponent::SteeringComponent(bool* iIsNetworkSending,
                                         int imagesPeriod,
                                         sem_t* varLock,
                                         ClientConnection* iClientConnection,
                                         vis::Control* iVisControl,
                                         lb::LBM* iLbm,
                                         net::Net * iNet,
                                         const topology::NetworkTopology *iNetTop,
                                         lb::SimulationState * iSimState) :
      net::PhasedBroadcast(iNet, iNetTop, iSimState, SPREADFACTOR), imagesPeriod(imagesPeriod),
          mLbm(iLbm), mSimState(iSimState), mVisControl(iVisControl)
    {
      Reset();
      AssignValues();
    }

    bool SteeringComponent::RequiresSeparateSteeringCore()
    {
      return false;
    }

    void SteeringComponent::ProgressFromChildren()
    {

    }
    void SteeringComponent::ProgressFromParent()
    {

    }
    void SteeringComponent::ProgressToChildren()
    {

    }
    void SteeringComponent::ProgressToParent()
    {

    }

    void SteeringComponent::PostReceiveFromChildren()
    {

    }
    void SteeringComponent::PostReceiveFromParent()
    {

    }

    void SteeringComponent::TopNodeAction()
    {

    }

    void SteeringComponent::Effect()
    {
      mSimState->DoRendering = ( (mSimState->TimeStep % imagesPeriod >= 0 && mSimState->TimeStep
          % imagesPeriod < (int) GetRoundTripLength()));
    }
  }
}
