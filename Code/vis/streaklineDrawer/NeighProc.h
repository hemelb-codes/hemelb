#ifndef HEMELB_VIS_NEIGHPROC_H
#define HEMELB_VIS_NEIGHPROC_H

#include <vector>

namespace hemelb
{
  namespace vis
  {
    namespace streaklinedrawer
    {
      struct NeighProc
      {
          proc_t id;
          site_t send_ps, recv_ps;
          site_t send_vs, recv_vs;

          std::vector<float> p_to_send, p_to_recv;
          float *v_to_send, *v_to_recv;

          site_t *s_to_send, *s_to_recv;
      };
    }
  }
}

#endif // HEMELB_VIS_NEIGHPROC_H
