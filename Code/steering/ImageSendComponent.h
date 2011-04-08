#ifndef HEMELB_STEERING_IMAGESENDCOMPONENT_H
#define HEMELB_STEERING_IMAGESENDCOMPONENT_H

#include "constants.h"
#include "lb/lb.h"
#include "lb/SimulationState.h"
#include "lb/LbmParameters.h"
#include "steering/ClientConnection.h"

namespace hemelb
{
  namespace steering
  {
    class ImageSendComponent
    {
      public:
        ImageSendComponent(lb::LBM* lbm,
                           lb::SimulationState* iSimState,
                           vis::Control* iControl,
                           const lb::LbmParameters* iLbmParams,
                           ClientConnection* iClientConnection);
        ~ImageSendComponent();

        void DoWork();

        bool ShouldRenderNewNetworkImage();

        bool isFrameReady;
        bool isConnected;
        int send_array_length;

      private:
        int SendSuccess(int iSocket, char * data, int length);

        ClientConnection* mClientConnection;
        lb::LBM* mLbm;
        lb::SimulationState* mSimState;
        vis::Control* mVisControl;
        const lb::LbmParameters* mLbmParams;

        char* xdrSendBuffer_pixel_data;
        double lastRender;

        // data per pixel inlude the data for the pixel location and 4 colours
        //in RGB format, thus #bytes per pixel are (sizeof(int)+4
        //rgb)=sizeof(int)+4*3*sizeof(unsigned char))
        static const int bytes_per_pixel_data = sizeof(int) + 4 * sizeof(unsigned char);
        // one int for colour_id and one for pixel id
        static const u_int pixel_data_bytes = vis::Screen::COLOURED_PIXELS_MAX
            * bytes_per_pixel_data;
        // it is assumed that the frame size is the only detail
        static const u_int frame_details_bytes = 1 * sizeof(int);
    };
  }
}

#endif /* HEMELB_STEERING_IMAGESENDCOMPONENT_H */
