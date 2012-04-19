#include "steering/ImageSendComponent.h"

namespace hemelb
{
  namespace steering
  {
    ImageSendComponent::ImageSendComponent(lb::SimulationState* iSimState,
                                           vis::Control* iControl,
                                           const lb::LbmParameters* iLbmParams,
                                           Network* iNetwork,
                                           unsigned int inletCountIn): inletCount(inletCountIn), MaxFramerate(25.0)
    {

    }

    ImageSendComponent::~ImageSendComponent()
    {
    }

    void ImageSendComponent::DoWork(const vis::PixelSet<vis::ResultPixel>* pix)
    {

    }

    bool ImageSendComponent::ShouldRenderNewNetworkImage()
    {
      return false;
    }
  }
}
