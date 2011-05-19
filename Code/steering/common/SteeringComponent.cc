#include "steering/SteeringComponent.h"

namespace hemelb
{
  namespace steering
  {
    void SteeringComponent::AssignValues()
    {
      mVisControl->mVisSettings.ctr_x += privateSteeringParams[0];
      mVisControl->mVisSettings.ctr_y += privateSteeringParams[1];
      mVisControl->mVisSettings.ctr_z += privateSteeringParams[2];

      float longitude = privateSteeringParams[3];
      float latitude = privateSteeringParams[4];

      float zoom = privateSteeringParams[5];

      mVisControl->mVisSettings.brightness = privateSteeringParams[6];

      // The minimum value here is by default 0.0 all the time
      mVisControl->mDomainStats.physical_velocity_threshold_max = privateSteeringParams[7];

      // The minimum value here is by default 0.0 all the time
      mVisControl->mDomainStats.physical_stress_threshold_max = privateSteeringParams[8];

      mVisControl->mDomainStats.physical_pressure_threshold_min = privateSteeringParams[9];
      mVisControl->mDomainStats.physical_pressure_threshold_max = privateSteeringParams[10];

      mVisControl->mVisSettings.glyphLength = privateSteeringParams[11];

      float pixels_x = privateSteeringParams[12];
      float pixels_y = privateSteeringParams[13];

      int newMouseX = int (privateSteeringParams[14]);
      int newMouseY = int (privateSteeringParams[15]);

      if (newMouseX != mVisControl->mVisSettings.mouse_x || newMouseY
          != mVisControl->mVisSettings.mouse_y)
      {
        updatedMouseCoords = true;
        mVisControl->mVisSettings.mouse_x = newMouseX;
        mVisControl->mVisSettings.mouse_y = newMouseY;
      }

      mSimState->SetIsTerminating(1 == (int) privateSteeringParams[16]);

      // To swap between glyphs and streak line rendering...
      // 0 - Only display the isosurfaces (wall pressure and stress)
      // 1 - Isosurface and glyphs
      // 2 - Wall pattern streak lines
      mVisControl->mVisSettings.mode = (vis::VisSettings::Mode) (privateSteeringParams[17]);

      mVisControl->mVisSettings.streaklines_per_pulsatile_period = privateSteeringParams[18];
      mVisControl->mVisSettings.streakline_length = privateSteeringParams[19];

      mSimState->SetDoRendering(1 == (int) privateSteeringParams[20]);

      mVisControl->UpdateImageSize((int) pixels_x, (int) pixels_y);

      distribn_t
          lattice_density_min =
              mLbm->ConvertPressureToLatticeUnits(mVisControl->mDomainStats.physical_pressure_threshold_min)
                  / Cs2;
      distribn_t
          lattice_density_max =
              mLbm->ConvertPressureToLatticeUnits(mVisControl->mDomainStats.physical_pressure_threshold_max)
                  / Cs2;
      distribn_t
          lattice_velocity_max =
              mLbm->ConvertVelocityToLatticeUnits(mVisControl->mDomainStats.physical_velocity_threshold_max);
      distribn_t
          lattice_stress_max =
              mLbm->ConvertStressToLatticeUnits(mVisControl->mDomainStats.physical_stress_threshold_max);

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

    void SteeringComponent::Reset()
    {
      isConnected = false;
      updatedMouseCoords = false;

      // scene center (dx,dy,dz)
      privateSteeringParams[0] = 0.0F;
      privateSteeringParams[1] = 0.0F;
      privateSteeringParams[2] = 0.0F;

      // longitude and latitude
      privateSteeringParams[3] = 45.0F;
      privateSteeringParams[4] = 45.0F;

      // zoom and brightness
      privateSteeringParams[5] = 1.0F;
      privateSteeringParams[6] = 0.03F;

      // velocity and stress ranges
      privateSteeringParams[7] = 0.1F;
      privateSteeringParams[8] = 0.1F;

      // Minimum pressure and maximum pressure for Colour mapping
      privateSteeringParams[9] = 80.0F;
      privateSteeringParams[10] = 120.0F;

      // Glyph length
      privateSteeringParams[11] = 1.0F;

      // Rendered frame size, pixel x and pixel y
      privateSteeringParams[12] = 512.0F;
      privateSteeringParams[13] = 512.0F;

      // x-y position of the mouse of the client
      privateSteeringParams[14] = -1.0F;
      privateSteeringParams[15] = -1.0F;

      // signal useful to terminate the simulation
      privateSteeringParams[16] = 0.0F;

      // Vis_mode
      privateSteeringParams[17] = 0.0F;

      // vis_streaklines_per_pulsatile_period
      privateSteeringParams[18] = 5.0F;

      // vis_streakline_length
      privateSteeringParams[19] = 100.0F;

      // Value of DoRendering
      privateSteeringParams[20] = 0.0F;
    }
  }
}
