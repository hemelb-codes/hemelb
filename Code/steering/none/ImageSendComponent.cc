#include "steering/ImageSendComponent.h"

namespace hemelb
{
  namespace steering
  {
    ImageSendComponent::ImageSendComponent(lb::LBM* lbm,
                                           lb::SimulationState* iSimState,
                                           vis::Control* iControl,
                                           const lb::LbmParameters* iLbmParams,
                                           Network* iNetwork)
    {

    }

    ImageSendComponent::~ImageSendComponent()
    {
    }

    void ImageSendComponent::DoWork(const vis::ScreenPixels* pix)
    {

    }

    bool ImageSendComponent::ShouldRenderNewNetworkImage()
    {
      return false;
    }
  }
}
