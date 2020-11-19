// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_STEERING_IMAGESENDCOMPONENT_H
#define HEMELB_STEERING_IMAGESENDCOMPONENT_H

#include <memory>
#include "constants.h"
#include "lb/SimulationState.h"
#include "lb/LbmParameters.h"
#include "steering/Network.h"
#include "steering/basic/SimulationParameters.h"
#include "vis/Control.h"

namespace hemelb
{
  namespace steering
  {
    class ImageSendComponent
    {
      public:
        ImageSendComponent(lb::SimulationState* iSimState, vis::Control* iControl,
                           const lb::LbmParameters* iLbmParams, Network* iNetwork,
                           unsigned inletCount);

        void DoWork(const vis::PixelSet<vis::ResultPixel>* pix);
        void SetMaxFramerate(float maxFramerate)
        {
          MaxFramerate = maxFramerate;
        }
        bool ShouldRenderNewNetworkImage();

        bool isConnected;
        int send_array_length;

      private:
        Network* mNetwork;
        lb::SimulationState* mSimState;
        vis::Control* mVisControl;
        const unsigned inletCount;
        float MaxFramerate;
        std::unique_ptr<char[]> xdrSendBuffer;
        double lastRender;

        // data per pixel is
        // 1 * int (pixel index)
        // 3 * int (pixel RGB)
        static const int bytes_per_pixel_data = 4 * 4;

        // Sent data:
        // 2 * int (pixelsX, pixelsY)
        // 1 * int (bytes of pixel data)
        // pixel data (variable, up to COLOURED_PIXELS_MAX * bytes_per_pixel_data)
        // SimulationParameters::paramsSizeB (metadata - mouse pressure and stress etc)
        static const unsigned int XdrIntLength = 4;
        static const unsigned int maxSendSize = 2 * XdrIntLength + 1 * XdrIntLength
            + vis::Screen::COLOURED_PIXELS_MAX * bytes_per_pixel_data
            + SimulationParameters::paramsSizeB;
    };
  }
}

#endif /* HEMELB_STEERING_IMAGESENDCOMPONENT_H */
