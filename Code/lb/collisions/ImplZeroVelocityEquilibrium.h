#ifndef HEMELB_LB_COLLISIONS_IMPLZEROVELOCITYEQUILIBRIUM_H
#define HEMELB_LB_COLLISIONS_IMPLZEROVELOCITYEQUILIBRIUM_H

#include "lb/collisions/WallCollision.h"

namespace hemelb
{
namespace lb
{
namespace collisions
{

class ImplZeroVelocityEquilibrium : public WallCollision
{
    void DoCollisions(double omega, int i, double *density, double *v_x, double *v_y,
      double *v_z, double f_neq[], Net* net);
};

}
}
}
#endif /* HEMELB_LB_COLLISIONS_IMPLZEROVELOCITYEQUILIBRIUM_H */
