// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "steering/ImageSendComponent.h"

namespace hemelb
{
  namespace steering
  {
    ImageSendComponent::ImageSendComponent(lb::SimulationState* iSimState, vis::Control* iControl,
                                           const lb::LbmParameters* iLbmParams, Network* iNetwork,
                                           unsigned int inletCountIn) :
        inletCount(inletCountIn), MaxFramerate(25.0)
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
