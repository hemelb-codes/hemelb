#include "steering/SteeringComponent.h"

namespace hemelb
{
  namespace steering
  {
    void SteeringComponent::AssignValues()
    {
      mVisControl->mVisSettings.ctr_x += privateSteeringParams[parameter::SceneCentreX];
      mVisControl->mVisSettings.ctr_y += privateSteeringParams[parameter::SceneCentreY];
      mVisControl->mVisSettings.ctr_z += privateSteeringParams[parameter::SceneCentreZ];

      float longitude = privateSteeringParams[parameter::Longitude];
      float latitude = privateSteeringParams[parameter::Latitude];

      float zoom = privateSteeringParams[parameter::Zoom];

      mVisControl->mVisSettings.brightness = privateSteeringParams[parameter::Brightness];

      // The minimum value here is by default 0.0 all the time
      mVisControl->mDomainStats.physical_velocity_threshold_max = privateSteeringParams[parameter::PhysicalVelocityThresholdMax];

      // The minimum value here is by default 0.0 all the time
      mVisControl->mDomainStats.physical_stress_threshold_max = privateSteeringParams[parameter::PhysicalStressThrehsholdMaximum];

      mVisControl->mDomainStats.physical_pressure_threshold_min = privateSteeringParams[parameter::PhysicalPressureThresholdMinimum];
      mVisControl->mDomainStats.physical_pressure_threshold_max = privateSteeringParams[10];

      mVisControl->mVisSettings.glyphLength = privateSteeringParams[parameter::GlyphLength];

      float pixels_x = privateSteeringParams[parameter::PixelsX];
      float pixels_y = privateSteeringParams[parameter::PixelsY];

      int newMouseX = int (privateSteeringParams[parameter::NewMouseX]);
      int newMouseY = int (privateSteeringParams[parameter::NewMouseY]);

      if (newMouseX != mVisControl->mVisSettings.mouse_x || newMouseY
          != mVisControl->mVisSettings.mouse_y)
      {
        updatedMouseCoords = true;
        mVisControl->mVisSettings.mouse_x = newMouseX;
        mVisControl->mVisSettings.mouse_y = newMouseY;
      }

      mSimState->SetIsTerminating(1 == (int) privateSteeringParams[parameter::SetIsTerminal]);

      // To swap between glyphs and streak line rendering...
      // 0 - Only display the isosurfaces (wall pressure and stress)
      // 1 - Isosurface and glyphs
      // 2 - Wall pattern streak lines
      mVisControl->mVisSettings.mode = (vis::VisSettings::Mode) (privateSteeringParams[parameter::Mode]);

      mVisControl->mVisSettings.streaklines_per_pulsatile_period = privateSteeringParams[parameter::StreaklinePerPulsatilePeriod];
      mVisControl->mVisSettings.streakline_length = privateSteeringParams[parameter::StreallineLength];

      mSimState->SetDoRendering(1 == (int) privateSteeringParams[parameter::SetDoRendering]);

      mVisControl->UpdateImageSize((int) pixels_x, (int) pixels_y);

      distribn_t
          lattice_density_min =
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

    void SteeringComponent::Reset(SimConfig& iSimConfig)
    {
      readyForNextImage = false;
      isConnected = false;
      updatedMouseCoords = false;

      // scene center (dx,dy,dz)
      privateSteeringParams[parameter::SceneCentreX] = 
	iSimConfig.VisCentre.x;
      privateSteeringParams[parameter::SceneCentreY] = 
	iSimConfig.VisCentre.y;
      privateSteeringParams[parameter::SceneCentreZ] = 
	iSimConfig.VisCentre.z;

      // longitude and latitude
      privateSteeringParams[parameter::Longitude] = 
	iSimConfig.VisLongitude;
      privateSteeringParams[parameter::Latitude] = 
	iSimConfig.VisLatitude;

      // zoom and brightness
      privateSteeringParams[parameter::Zoom] = 
	iSimConfig.VisZoom;
      privateSteeringParams[parameter::Brightness] = 
	iSimConfig.VisBrightness;

      // velocity and stress ranges
      privateSteeringParams[parameter::PhysicalVelocityThresholdMax] = iSimConfig.MaxVelocity;
      privateSteeringParams[parameter::PhysicalStressThrehsholdMaximum] = iSimConfig.MaxStress;

      // Minimum pressure and maximum pressure for Colour mapping
      privateSteeringParams[parameter::PhysicalPressureThresholdMinimum] = 80.0F;
      privateSteeringParams[10] = 120.0F;

      // Glyph length
      privateSteeringParams[parameter::GlyphLength] = 1.0F;

      // Rendered frame size, pixel x and pixel y
      privateSteeringParams[parameter::PixelsX] = 512.0F;
      privateSteeringParams[parameter::PixelsY] = 512.0F;

      // x-y position of the mouse of the client
      privateSteeringParams[parameter::NewMouseX] = -1.0F;
      privateSteeringParams[parameter::NewMouseY] = -1.0F;

      // signal useful to terminate the simulation
      privateSteeringParams[parameter::SetIsTerminal] = 0.0F;

      // Vis_mode
      privateSteeringParams[parameter::Mode] = 0.0F;

      // vis_streaklines_per_pulsatile_period
      privateSteeringParams[parameter::StreaklinePerPulsatilePeriod] = 5.0F;

      // vis_streakline_length
      privateSteeringParams[parameter::StreallineLength] = 100.0F;

      // Value of DoRendering
      privateSteeringParams[parameter::SetDoRendering] = 0.0F;
    }
  }
}
