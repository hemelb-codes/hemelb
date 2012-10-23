// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_VIS_STREAKLINEDRAWER_STREAKPIXEL_H
#define HEMELB_VIS_STREAKLINEDRAWER_STREAKPIXEL_H

#include "log/Logger.h"
#include "vis/BasicPixel.h"

namespace hemelb
{
  namespace vis
  {
    namespace streaklinedrawer
    {
      class StreakPixel : public BasicPixel
      {
        public:
          StreakPixel()
          {

          }

          StreakPixel(int i, int j, float particleVelocity, float particleZ, int particleInlet) :
            BasicPixel(i, j), particle_vel(particleVelocity), particle_z(particleZ), particle_inlet_id(particleInlet)
          {

          }

          void Combine(const StreakPixel& inPixel)
          {
            if (inPixel.particle_z < particle_z)
            {
              particle_z = inPixel.particle_z;
              particle_vel = inPixel.particle_vel;
              particle_inlet_id = inPixel.particle_inlet_id;
            }
          }

          float GetParticleVelocity() const
          {
            return particle_vel;
          }

          /**
           * Produces an MPI Datatype object but doesn't commit it or manage its memory.
           * @return
           */
          static MPI_Datatype GetMPIType()
          {
            const int typeCount = 5;
            int blocklengths[typeCount] = { 1, 1, 1, 1, 1 };

            MPI_Datatype types[typeCount] = { MPI_INT, MPI_INT, MPI_FLOAT, MPI_FLOAT, MPI_INT };

            StreakPixel example;

            MPI_Aint displacements[typeCount];

            MPI_Get_address(&example.i, &displacements[0]);
            MPI_Get_address(&example.j, &displacements[1]);
            MPI_Get_address(&example.particle_vel, &displacements[2]);
            MPI_Get_address(&example.particle_z, &displacements[3]);
            MPI_Get_address(&example.particle_inlet_id, &displacements[4]);

            for (int ii = typeCount - 1; ii >= 0; --ii)
            {
              displacements[ii] -= displacements[0];
            }

            MPI_Datatype ret;

            MPI_Type_struct(typeCount, blocklengths, displacements, types, &ret);

            return ret;
          }

          void LogDebuggingInformation() const
          {
            log::Logger::Log<log::Trace, log::OnePerCore>("Streak pixel at (%i,%i) with "
              "(source inlet, velocity, z) = (%d, %f, %f)", GetI(), GetJ(), particle_inlet_id, particle_vel, particle_z);
          }

        private:
          float particle_vel;
          float particle_z;
          int particle_inlet_id;
      };
    }
  }
}

#endif /* HEMELB_VIS_STREAKLINEDRAWER_STREAKPIXEL_H */
