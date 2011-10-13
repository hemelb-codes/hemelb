#ifndef HEMELB_VIS_RAYTRACER_RAYTRACER_RAYPIXEL_H
#define HEMELB_VIS_RAYTRACER_RAYTRACER_RAYPIXEL_H

#include "util/utilityFunctions.h"
#include "vis/BasicPixel.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer
    {
      class RayPixel : public BasicPixel
      {
        public:
          RayPixel()
          {

          }

          RayPixel(int i,
                   int j,
                   float t,
                   float dt,
                   float density,
                   float stress,
                   float velArrIn[3],
                   float stressArrIn[3]) :
            BasicPixel(i, j), t(t), dt(dt), density(density), stress(stress)
          {
            for (int ii = 0; ii < 3; ++ii)
            {
              velArr[ii] = velArrIn[ii];
              stressArr[ii] = stressArrIn[ii];
            }
          }

          void Combine(const RayPixel& other)
          {
            // Both are ray-tracing
            for (int ii = 0; ii < 3; ++ii)
            {
              velArr[ii] += other.velArr[ii];
              stressArr[ii] += other.stressArr[ii];
            }

            dt += other.dt;

            if (other.t < t)
            {
              t = other.t;
              density = other.density;
              stress = other.stress;
            }
          }

          float GetT() const
          {
            return t;
          }

          float GetDT() const
          {
            return dt;
          }

          float GetDensity() const
          {
            return density;
          }

          float GetStress() const
          {
            return stress;
          }

          const float* GetVelArray() const
          {
            return velArr;
          }

          const float* GetStressArray() const
          {
            return stressArr;
          }

          static void PickColour(float value, float colour[3])
          {
            colour[0] = util::NumericalFunctions::enforceBounds<float>(4.F * value - 2.F, 0.F, 1.F);
            colour[1] = util::NumericalFunctions::enforceBounds<float>(2.F - 4.F
                * (float) fabs(value - 0.5F), 0.F, 1.F);
            colour[2] = util::NumericalFunctions::enforceBounds<float>(2.F - 4.F * value, 0.F, 1.F);
          }

          /**
           * Produces an MPI Datatype object but doesn't commit it or manage its memory.
           * @return
           */
          static MPI_Datatype GetMPIType()
          {
            const int typeCount = 8;
            int blocklengths[typeCount];
            MPI_Datatype types[typeCount];

            blocklengths[0] = blocklengths[1] = 1;
            types[0] = types[1] = MpiDataType<int> ();

            for (int ii = 2; ii < typeCount; ++ii)
            {
              blocklengths[ii] = ii < 6
                ? 1
                : 3;
              types[ii] = MPI_FLOAT;
            }

            MPI_Aint displacements[typeCount];
            RayPixel example;

            MPI_Get_address(&example.i, &displacements[0]);
            MPI_Get_address(&example.j, &displacements[1]);
            MPI_Get_address(&example.t, &displacements[2]);
            MPI_Get_address(&example.dt, &displacements[3]);
            MPI_Get_address(&example.density, &displacements[4]);
            MPI_Get_address(&example.stress, &displacements[5]);
            MPI_Get_address(&example.velArr, &displacements[6]);
            MPI_Get_address(&example.stressArr, &displacements[7]);

            for (int ii = typeCount - 1; ii >= 0; --ii)
            {
              MPI_Aint baseDisplacement;
              MPI_Get_address(&example, &baseDisplacement);
              displacements[ii] -= baseDisplacement;
            }

            MPI_Datatype ret;

            MPI_Type_struct(typeCount, blocklengths, displacements, types, &ret);

            return ret;
          }

        private:
          float t, dt;
          float density;
          float stress;
          float velArr[3];
          float stressArr[3];
      };
    }
  }
}

#endif //HEMELB_VIS_RAYTRACER_RAYTRACER_RAYPIXEL_H
