#ifndef HEMELB_VIS_DOMAINSTATS_H
#define HEMELB_VIS_DOMAINSTATS_H

namespace hemelb
{
  namespace vis
  {
    struct DomainStats
    {
      public:
        float velocity_threshold_max_inv;
        float stress_threshold_max_inv;
        float density_threshold_min, density_threshold_minmax_inv;
        float physical_pressure_threshold_min;
        float physical_pressure_threshold_max;
        float physical_velocity_threshold_max;
        float physical_stress_threshold_max;
    };
  }
}

#endif /* HEMELB_VIS_DOMAINSTATS_H */
