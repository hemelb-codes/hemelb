#include "steering/ImageSendComponent.h"

namespace hemelb
{
  namespace steering
  {
    ImageSendComponent::ImageSendComponent(lb::LBM* lbm,
                                           lb::SimulationState* iSimState,
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

    double ImageSendComponent::frameTiming(void)
    {

    }

    int ImageSendComponent::SendSuccess(int iSocket, char * data, int length)
    {

    }

    bool ImageSendComponent::ShouldRenderNewNetworkImage()
    {
      return false;
    }
  }
}
