#ifndef HEMELB_LB_COLLISIONS_IMPLZEROVELOCITYBOUNDARYDENSITY_H
#define HEMELB_LB_COLLISIONS_IMPLZEROVELOCITYBOUNDARYDENSITY_H

#include "lb/collisions/InletOutletWallCollision.h"

namespace hemelb
{
namespace lb
{
namespace collisions
{

class ImplZeroVelocityBoundaryDensity : public InletOutletWallCollision
{
  public:
    ImplZeroVelocityBoundaryDensity(double* iBoundaryDensityArray);
    void DoCollisions(double omega, int i, double *density, double *v_x, double *v_y,
                      double *v_z, double f_neq[], Net* net);

  private:
    double* mBoundaryDensityArray;
};

}
}
}

#endif /* HEMELB_LB_COLLISIONS_IMPLZEROVELOCITYBOUNDARYDENSITY_H */
