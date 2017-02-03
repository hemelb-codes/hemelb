// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include "steering/SteeringComponent.h"

namespace hemelb
{
  namespace steering
  {
    void SteeringComponent::AssignValues()
    {
      mVisControl->visSettings.ctr_x += privateSteeringParams[SceneCentreX];
      mVisControl->visSettings.ctr_y += privateSteeringParams[SceneCentreY];
      mVisControl->visSettings.ctr_z += privateSteeringParams[SceneCentreZ];

      float longitude = privateSteeringParams[Longitude];
      float latitude = privateSteeringParams[Latitude];

      float zoom = privateSteeringParams[Zoom];

      mVisControl->visSettings.brightness = privateSteeringParams[Brightness];

      // The minimum value here is by default 0.0 all the time
      mVisControl->domainStats.physical_velocity_threshold_max =
          privateSteeringParams[PhysicalVelocityThresholdMax];

      // The minimum value here is by default 0.0 all the time
      mVisControl->domainStats.physical_stress_threshold_max =
          privateSteeringParams[PhysicalStressThrehsholdMaximum];

      mVisControl->domainStats.physical_pressure_threshold_min =
          privateSteeringParams[PhysicalPressureThresholdMinimum];
      mVisControl->domainStats.physical_pressure_threshold_max = privateSteeringParams[10];

      mVisControl->visSettings.glyphLength = privateSteeringParams[GlyphLength];

      float pixels_x = privateSteeringParams[PixelsX];
      float pixels_y = privateSteeringParams[PixelsY];

      int newMouseX = int(privateSteeringParams[NewMouseX]);
      int newMouseY = int(privateSteeringParams[NewMouseY]);

      if (newMouseX != mVisControl->visSettings.mouse_x
          || newMouseY != mVisControl->visSettings.mouse_y)
      {
        updatedMouseCoords = true;
        mVisControl->visSettings.mouse_x = newMouseX;
        mVisControl->visSettings.mouse_y = newMouseY;
      }

      mSimState->SetIsTerminating(1 == (int) privateSteeringParams[SetIsTerminal]);

      // To swap between glyphs and streak line rendering...
      // 0 - Only display the isosurfaces (wall pressure and stress)
      // 1 - Isosurface and glyphs
      // 2 - Wall pattern streak lines
      mVisControl->visSettings.mode = (vis::VisSettings::Mode) (privateSteeringParams[Mode]);

      mVisControl->visSettings.streaklines_per_simulation =
          privateSteeringParams[StreaklinePerSimulation];
      mVisControl->visSettings.streakline_length = privateSteeringParams[StreaklineLength];

      mSimState->SetIsRendering(1 == (int) privateSteeringParams[SetDoRendering]);
      if (mSimState->IsRendering())
      {
        readyForNextImage = false;
      }

      mVisControl->UpdateImageSize((int) pixels_x, (int) pixels_y);

      distribn_t lattice_density_min =
          mUnits->ConvertPressureToLatticeUnits(mVisControl->domainStats.physical_pressure_threshold_min)
              / Cs2;
      distribn_t lattice_density_max =
          mUnits->ConvertPressureToLatticeUnits(mVisControl->domainStats.physical_pressure_threshold_max)
              / Cs2;
      distribn_t lattice_velocity_max =
          mUnits->ConvertVelocityToLatticeUnits(mVisControl->domainStats.physical_velocity_threshold_max);
      distribn_t lattice_stress_max =
          mUnits->ConvertStressToLatticeUnits(mVisControl->domainStats.physical_stress_threshold_max);

      mVisControl->SetProjection((int) pixels_x,
                                 (int) pixels_y,
                                 mVisControl->visSettings.ctr_x,
                                 mVisControl->visSettings.ctr_y,
                                 mVisControl->visSettings.ctr_z,
                                 longitude,
                                 latitude,
                                 zoom);
      if (imageSendComponent != nullptr)
      {
        imageSendComponent->SetMaxFramerate(privateSteeringParams[MaxFramerate]);
      }
      mVisControl->domainStats.density_threshold_min = lattice_density_min;
      mVisControl->domainStats.density_threshold_minmax_inv = 1.0F
          / (lattice_density_max - lattice_density_min);
      mVisControl->domainStats.velocity_threshold_max_inv = 1.0F / lattice_velocity_max;
      mVisControl->domainStats.stress_threshold_max_inv = 1.0F / lattice_stress_max;
    }

    void SteeringComponent::ClearValues()
    {
      readyForNextImage = false;
      isConnected = false;
      updatedMouseCoords = false;

      // scene center (dx,dy,dz)
      privateSteeringParams[SceneCentreX] = simConfig->GetVisualisationCentre().x;
      privateSteeringParams[SceneCentreY] = simConfig->GetVisualisationCentre().y;
      privateSteeringParams[SceneCentreZ] = simConfig->GetVisualisationCentre().z;

      // longitude and latitude
      privateSteeringParams[Longitude] = simConfig->GetVisualisationLongitude();
      privateSteeringParams[Latitude] = simConfig->GetVisualisationLatitude();

      // zoom and brightness
      privateSteeringParams[Zoom] = simConfig->GetVisualisationZoom();
      privateSteeringParams[Brightness] = simConfig->GetVisualisationBrightness();

      // velocity and stress ranges
      privateSteeringParams[PhysicalVelocityThresholdMax] = simConfig->GetMaximumVelocity();
      privateSteeringParams[PhysicalStressThrehsholdMaximum] = simConfig->GetMaximumStress();

      // Minimum pressure and maximum pressure for Colour mapping
      privateSteeringParams[PhysicalPressureThresholdMinimum] = 80.0F;
      privateSteeringParams[10] = 120.0F;

      // Glyph length
      privateSteeringParams[GlyphLength] = 1.0F;

      // Rendered frame size, pixel x and pixel y
      privateSteeringParams[PixelsX] = 512.0F;
      privateSteeringParams[PixelsY] = 512.0F;

      // x-y position of the mouse of the client
      privateSteeringParams[NewMouseX] = -1.0F;
      privateSteeringParams[NewMouseY] = -1.0F;

      // signal useful to terminate the simulation
      privateSteeringParams[SetIsTerminal] = 0.0F;

      // Vis_mode
      privateSteeringParams[Mode] = 0.0F;

      // vis_streaklines_per_pulsatile_period
      privateSteeringParams[StreaklinePerSimulation] = 5.0F;

      // vis_streakline_length
      privateSteeringParams[StreaklineLength] = 100.0F;

      privateSteeringParams[MaxFramerate] = 25.0F;

      // Value of DoRendering
      privateSteeringParams[SetDoRendering] = 0.0F;
    }
  }
}
