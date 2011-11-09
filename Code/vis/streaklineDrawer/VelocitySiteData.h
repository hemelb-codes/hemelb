#ifndef HEMELB_VIS_STREAKLINEDRAWER_VELOCITYSITEDATA_H
#define HEMELB_VIS_STREAKLINEDRAWER_VELOCITYSITEDATA_H

namespace hemelb
{
  namespace vis
  {
    namespace streaklinedrawer
    {
      //Class for information about the velocity field at some point.
      class VelocitySiteData
      {
        public:
          VelocitySiteData()
          {
            proc_id = -1;
            counter = 0;
          }

          proc_t proc_id;
          site_t counter, site_id;
          float vx, vy, vz;
      };
    }
  }
}

#endif //HEMELB_VIS_STREAKLINEDRAWER_VELOCITYSITEDATA_H
