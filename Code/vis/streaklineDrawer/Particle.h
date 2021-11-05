// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_VIS_STREAKLINEDRAWER_PARTICLE_H
#define HEMELB_VIS_STREAKLINEDRAWER_PARTICLE_H

#include "util/Vector3D.h"

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

          util::Vector3D<float> position;
          util::Vector3D<float> velocity;
          float vel;
          unsigned int inletID;
      };

    }
  }
  namespace net
  {
    template<>
    MPI_Datatype MpiDataTypeTraits<hemelb::vis::streaklinedrawer::Particle>::RegisterMpiDataType();
  }
}

#endif // HEMELB_VIS_STREAKLINEDRAWER_PARTICLE_H
