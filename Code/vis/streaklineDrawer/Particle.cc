// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include "net/mpi.h"
#include "vis/streaklineDrawer/Particle.h"
#include "constants.h"

namespace hemelb
{
  namespace vis
  {
    namespace streaklinedrawer
    {
      Particle::Particle() :
          position(NO_VALUE), velocity(NO_VALUE), vel(NO_VALUE), inletID()
      {

      }

      Particle::Particle(float iX, float iY, float iZ, unsigned int iInletId) :
          position(iX, iY, iZ), velocity(0), vel(0), inletID(iInletId)
      {
      }
    }
  }

  namespace net
  {
    template<>
    MPI_Datatype MpiDataTypeTraits<hemelb::vis::streaklinedrawer::Particle>::RegisterMpiDataType()
    {
      HEMELB_MPI_TYPE_BEGIN(type, vis::streaklinedrawer::Particle, 3);

      HEMELB_MPI_TYPE_ADD_MEMBER(position);
      HEMELB_MPI_TYPE_ADD_MEMBER(vel);
      HEMELB_MPI_TYPE_ADD_MEMBER(inletID);

      HEMELB_MPI_TYPE_END(type, vis::streaklinedrawer::Particle);

      HEMELB_MPI_CALL(MPI_Type_commit, (&type));
      return type;
    }
  }
}
