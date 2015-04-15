// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_STEERING_IMAGESENDCOMPONENT_H
#define HEMELB_STEERING_IMAGESENDCOMPONENT_H

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
        ~ImageSendComponent();

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
        char* xdrSendBuffer;
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
