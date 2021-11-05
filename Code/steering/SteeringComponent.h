// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_STEERING_STEERINGCOMPONENT_H
#define HEMELB_STEERING_STEERINGCOMPONENT_H

#include "net/PhasedBroadcastRegular.h"
#include "lb/SimulationState.h"
#include "steering/Network.h"
#include "vis/DomainStats.h"

namespace hemelb
{
  namespace configuration {
    class SimConfig;
  }
  namespace util {
    class UnitConverter;
  }
  namespace vis {
    class Control;
  }

  namespace steering
  {
    class ImageSendComponent;

    enum parameter
    {
      SceneCentreX = 0,
      SceneCentreY = 1,
      SceneCentreZ = 2,
      Longitude = 3,
      Latitude = 4,
      Zoom = 5,
      Brightness = 6,
      PhysicalVelocityThresholdMax = 7,
      PhysicalStressThrehsholdMaximum = 8,
      PhysicalPressureThresholdMinimum = 9,
      PhysicalPressureThresholdMaximum = 10,
      GlyphLength = 11,
      PixelsX = 12,
      PixelsY = 13,
      NewMouseX = 14,
      NewMouseY = 15,
      SetIsTerminal = 16,
      Mode = 17,
      StreaklinePerSimulation = 18,
      StreaklineLength = 19,
      MaxFramerate = 20,
      SetDoRendering = 21
    };

    /**
     * SteeringComponent - class for passing steering data to all nodes.
     *
     * We pass this data at regular intervals. No initial action is required by all nodes, and
     * we only need to pass from the top-most node (which handles network communication) downwards,
     * on one iteration between each pair of consecutive depths.
     */
    class SteeringComponent : public net::PhasedBroadcastRegular<false, 1, 0, true, false>
    {
      public:
        SteeringComponent(Network* iNetwork, vis::Control* iVisControl,
                          steering::ImageSendComponent* imageSendComponent, net::Net * iNet,
                          lb::SimulationState * iSimState, configuration::SimConfig* iSimConfig,
                          const util::UnitConverter* iUnits);

        static bool RequiresSeparateSteeringCore();

        /*
         * This function initialises all of the steering parameters, on all nodes.
         * Although this appears to be part of the removed post-unstable reset functionality, it is also used during initialisation
         * so is retained #244
         */
        void ClearValues();

        bool readyForNextImage;
        bool updatedMouseCoords;

      protected:
        void ProgressFromParent(unsigned long splayNumber);
        void ProgressToChildren(unsigned long splayNumber);

        void TopNodeAction();
        void Effect();

      private:
        void AssignValues();

        const static int STEERABLE_PARAMETERS = 21;
        const static unsigned int SPREADFACTOR = 10;

        bool isConnected;

        Network* mNetwork;
        lb::SimulationState* mSimState;
        vis::Control* mVisControl;
        steering::ImageSendComponent* imageSendComponent;
        float privateSteeringParams[STEERABLE_PARAMETERS + 1];
        const util::UnitConverter* mUnits;
        configuration::SimConfig* simConfig;
    };
  }
}

#endif /* HEMELB_STEERING_STEERINGCOMPONENT_H */
