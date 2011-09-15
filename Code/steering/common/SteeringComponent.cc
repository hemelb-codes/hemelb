#include "steering/SteeringComponent.h"

namespace hemelb
{
  namespace steering
  {
    void SteeringComponent::AssignValues()
    {
      mVisControl->mVisSettings.ctr_x += privateSteeringParams[SceneCentreX];
      mVisControl->mVisSettings.ctr_y += privateSteeringParams[SceneCentreY];
      mVisControl->mVisSettings.ctr_z += privateSteeringParams[SceneCentreZ];

      float longitude = privateSteeringParams[Longitude];
      float latitude = privateSteeringParams[Latitude];

      float zoom = privateSteeringParams[Zoom];

      mVisControl->mVisSettings.brightness = privateSteeringParams[Brightness];

      // The minimum value here is by default 0.0 all the time
      mVisControl->mDomainStats.physical_velocity_threshold_max
          = privateSteeringParams[PhysicalVelocityThresholdMax];

      // The minimum value here is by default 0.0 all the time
      mVisControl->mDomainStats.physical_stress_threshold_max
          = privateSteeringParams[PhysicalStressThrehsholdMaximum];

      mVisControl->mDomainStats.physical_pressure_threshold_min
          = privateSteeringParams[PhysicalPressureThresholdMinimum];
      mVisControl->mDomainStats.physical_pressure_threshold_max = privateSteeringParams[10];

      mVisControl->mVisSettings.glyphLength = privateSteeringParams[GlyphLength];

      float pixels_x = privateSteeringParams[PixelsX];
      float pixels_y = privateSteeringParams[PixelsY];

      int newMouseX = int (privateSteeringParams[NewMouseX]);
      int newMouseY = int (privateSteeringParams[NewMouseY]);

      if (newMouseX != mVisControl->mVisSettings.mouse_x || newMouseY
          != mVisControl->mVisSettings.mouse_y)
      {
        updatedMouseCoords = true;
        mVisControl->mVisSettings.mouse_x = newMouseX;
        mVisControl->mVisSettings.mouse_y = newMouseY;
      }

      mSimState->SetIsTerminating(1 == (int) privateSteeringParams[SetIsTerminal]);

      // To swap between glyphs and streak line rendering...
      // 0 - Only display the isosurfaces (wall pressure and stress)
      // 1 - Isosurface and glyphs
      // 2 - Wall pattern streak lines
      mVisControl->mVisSettings.mode
          = (vis::VisSettings::Mode) (privateSteeringParams[Mode]);

      mVisControl->mVisSettings.streaklines_per_pulsatile_period
          = privateSteeringParams[StreaklinePerPulsatilePeriod];
      mVisControl->mVisSettings.streakline_length
          = privateSteeringParams[StreallineLength];

      mSimState->SetDoRendering(1 == (int) privateSteeringParams[SetDoRendering]);

      mVisControl->UpdateImageSize((int) pixels_x, (int) pixels_y);

      distribn_t lattice_density_min =
          mUnits->ConvertPressureToLatticeUnits(mVisControl->mDomainStats.physical_pressure_threshold_min)
                  / Cs2;
      distribn_t
          lattice_density_max =
              mUnits->ConvertPressureToLatticeUnits(mVisControl->mDomainStats.physical_pressure_threshold_max)
                  / Cs2;
      distribn_t
          lattice_velocity_max =
              mUnits->ConvertVelocityToLatticeUnits(mVisControl->mDomainStats.physical_velocity_threshold_max);
      distribn_t
          lattice_stress_max =
              mUnits->ConvertStressToLatticeUnits(mVisControl->mDomainStats.physical_stress_threshold_max);

      mVisControl->SetProjection((int) pixels_x,
                                 (int) pixels_y,
                                 mVisControl->mVisSettings.ctr_x,
                                 mVisControl->mVisSettings.ctr_y,
                                 mVisControl->mVisSettings.ctr_z,
                                 longitude,
                                 latitude,
                                 zoom);

      mVisControl->mDomainStats.density_threshold_min = lattice_density_min;
      mVisControl->mDomainStats.density_threshold_minmax_inv = 1.0F / (lattice_density_max
          - lattice_density_min);
      mVisControl->mDomainStats.velocity_threshold_max_inv = 1.0F / lattice_velocity_max;
      mVisControl->mDomainStats.stress_threshold_max_inv = 1.0F / lattice_stress_max;
    }

    void SteeringComponent::Reset(SimConfig* iSimConfig)
    {
      readyForNextImage = false;
      isConnected = false;
      updatedMouseCoords = false;

      // scene center (dx,dy,dz)
      privateSteeringParams[SceneCentreX] = iSimConfig->VisCentre.x;
      privateSteeringParams[SceneCentreY] = iSimConfig->VisCentre.y;
      privateSteeringParams[SceneCentreZ] = iSimConfig->VisCentre.z;

      // longitude and latitude
      privateSteeringParams[Longitude] = iSimConfig->VisLongitude;
      privateSteeringParams[Latitude] = iSimConfig->VisLatitude;

      // zoom and brightness
      privateSteeringParams[Zoom] = iSimConfig->VisZoom;
      privateSteeringParams[Brightness] = iSimConfig->VisBrightness;

      // velocity and stress ranges
      privateSteeringParams[PhysicalVelocityThresholdMax] = iSimConfig->MaxVelocity;
      privateSteeringParams[PhysicalStressThrehsholdMaximum] = iSimConfig->MaxStress;

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
      privateSteeringParams[StreaklinePerPulsatilePeriod] = 5.0F;

      // vis_streakline_length
      privateSteeringParams[StreallineLength] = 100.0F;

      // Value of DoRendering
      privateSteeringParams[SetDoRendering] = 0.0F;
    }
  }
}
