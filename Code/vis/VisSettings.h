#ifndef HEMELB_VIS_VISSETTINGS_H
#define HEMELB_VIS_VISSETTINGS_H

#include "lb/LbmParameters.h"

namespace hemelb
{
  namespace vis
  {
    struct VisSettings
    {
        // better public member vars than globals!
        int mode;

        float ctr_x, ctr_y, ctr_z;
        float streaklines_per_pulsatile_period, streakline_length;
        double mouse_pressure, mouse_stress;
        float brightness;

        lb::StressTypes mStressType;

        int mouse_x, mouse_y;
    };
  }
}

#endif /* HEMELB_VIS_VISSETTINGS_H */
