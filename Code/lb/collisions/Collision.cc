#include "lb/collisions/Collision.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      bool Collision::PostStep(double omega, int i, double *density, double *v_x,
        double *v_y, double *v_z, double f_neq[], Net* net)
      {
        // Standard implementation - do nothing. Return false to show that we've
        // done nothing.
        return false;
      }

      unsigned int Collision::GetBoundaryConfig(Net* net, int i)
      {
        return (net->net_site_data[i] & BOUNDARY_CONFIG_MASK) >> BOUNDARY_CONFIG_SHIFT;
      }
    }
  }
}
