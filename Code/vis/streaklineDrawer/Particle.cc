#include "Particle.h"

namespace hemelb
{
  namespace vis
  {
    namespace streaklinedrawer
    {
      Particle::Particle()
      {

      }

      Particle::Particle(float iX, float iY, float iZ, unsigned int iInletId) :
        x(iX), y(iY), z(iZ), vx(0), vy(0), vz(0), vel(0), inletID(iInletId)
      {
      }

      Particle::Particle(float iX, float iY, float iZ, float iVel, unsigned int iInletId) :
        x(iX), y(iY), z(iZ), vx(0), vy(0), vz(0), vel(iVel), inletID(iInletId)
      {
      }

    }
  }
}
