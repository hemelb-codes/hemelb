#include "steering/ImageSendComponent.h"

namespace hemelb
{
  namespace steering
  {
    ImageSendComponent::ImageSendComponent(lb::LBM* lbm,
                                           lb::SimulationState* iSimState,
                                           vis::Control* iControl,
                                           const lb::LbmParameters* iLbmParams,
                                           ClientConnection* iClientConnection)
    {

    }

    ImageSendComponent::~ImageSendComponent()
    {
    }

    void ImageSendComponent::DoWork()
    {

    }

    int ImageSendComponent::SendSuccess(int iSocket, char * data, int length)
    {
      return -1;
    }

    bool ImageSendComponent::ShouldRenderNewNetworkImage()
    {
      return false;
    }
  }
}
