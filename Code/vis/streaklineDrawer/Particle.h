#ifndef HEMELB_VIS_PARTICLE_H
#define HEMELB_VIS_PARTICLE_H

namespace hemelb
{
  namespace vis
  {
    namespace streaklinedrawer
    {
      struct Particle
      {
        public:
          Particle();

          Particle(float iX, float iY, float iZ, unsigned int iInletId);

          Particle(float iX, float iY, float iZ, float iVel, unsigned int iInletId);

          float x, y, z;
          float vx, vy, vz;
          float vel;
          unsigned int inletID;
      };

    }
  }
}

#endif // HEMELB_VIS_PARTICLE_H
